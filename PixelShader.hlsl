#define RANDOM_IA 16807
#define RANDOM_IM 2147483647
#define RANDOM_AM (1.0f/float(RANDOM_IM))
#define RANDOM_IQ 127773u
#define RANDOM_IR 2836
#define RANDOM_MASK 123459876
#define MAX_UINT 4294967290u
#define MAX_INT 2147483647

static float3 white = float3(1.0, 1.0, 1.0);
static float3 black = float3(0.0, 0.0, 0.0);
static float3 light_blue = float3(0.5, 0.7, 1.0);
struct VS_OUTPUT
{
	float4 pos: SV_POSITION;
	float4 color: COLOR;
	float4 p : POSITION;
};

cbuffer ConstantBuffer : register(b0)
{
	float3 origin;
	float width;
	float3 horizontal;
	float height;
	float3 vertical;
	int numSpheres;
	float3 center;
	int numBounces;
	int numSamples;
	int randomVal;
	float aspectRatio;
};
struct NumberGenerator {
	int seed; // Used to generate values.

			  // Returns the current random float.
	float GetCurrentFloat() {
		Cycle();
		return RANDOM_AM * seed;
	}

	// Returns the current random int.
	int GetCurrentInt() {
		Cycle();
		return seed;
	}

	// Generates the next number in the sequence.
	void Cycle() {
		seed ^= RANDOM_MASK;
		int k = seed / RANDOM_IQ;
		seed = RANDOM_IA * (seed - k * RANDOM_IQ) - RANDOM_IR * k;

		if (seed < 0)
			seed += RANDOM_IM;

		seed ^= RANDOM_MASK;
	}

	// Cycles the generator based on the input count. Useful for generating a thread unique seed.
	// PERFORMANCE - O(N)
	void Cycle(const uint _count) {
		for (uint i = 0; i < _count; ++i)
			Cycle();
	}

	// Returns a random float within the input range.
	float GetRandomFloat(const float low, const float high) {
		float v = GetCurrentFloat();
		return low * (1.0f - v) + high * v;
	}

	// Sets the seed
	void SetSeed(const uint value) {
		seed = int(value);
		Cycle();
	}
};
static NumberGenerator rd;

float random(float low, float high)
{
	return rd.GetRandomFloat(low, high);
}

float rand_1_05(in float2 uv)
{
	float2 noise = (frac(sin(dot(uv, float2(12.9898, 78.233)*2.0)) * 43758.5453));
	return abs(noise.x + noise.y) * 0.5;
}

float2 rand_2_10(in float2 uv) {
	float noiseX = (frac(sin(dot(uv, float2(12.9898, 78.233) * 2.0)) * 43758.5453));
	float noiseY = sqrt(1 - noiseX * noiseX);
	return float2(noiseX, noiseY);
}
float2 rand_2_0004(in float2 uv)
{
	float noiseX = (frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453));
	float noiseY = (frac(sin(dot(uv, float2(12.9898, 78.233) * 2.0)) * 43758.5453));
	return float2(noiseX, noiseY) * 0.004;
}

uint hash(uint x) {
	x += (x << 10u);
	x ^= (x >> 6u);
	x += (x << 3u);
	x ^= (x >> 11u);
	x += (x << 15u);
	return x;
}
float floatConstruct(uint m) {
	const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
	const uint ieeeOne = 0x3F800000u; // 1.0 in IEEE binary32

	m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
	m |= ieeeOne;                          // Add fractional part to 1.0

	float  f = asfloat(m);                 // Range [1:2]
	return f - 1.0;                        // Range [0:1]
}

uint hash(uint3 v) { return hash(v.x ^ hash(v.y) ^ hash(v.z)); }

float rndSd(float3  v) { return floatConstruct(hash(asuint(v))); }

float3 random_in_unit_sphere()
{
	float3 result;
	do {
		result = 2 * float3(random(0.0, 1.0),random(0.0,1.0),random(0.0,1.0)) - float3(1.0, 1.0, 1.0);
	} while (dot(result, result) > 1.0);
	return result;
}

struct Ray
{
	float3 origin;
	float3 direction;
	float3 time;

	float3 point_at_parameter(float t)
	{
		return origin + t*direction;
	}
};

struct Sphere
{
	float3 startPos;
	float radius;
	float3 material;
	float matAdd;
	float3 endPos;
	int matType;
	float3 origin(float t)
	{
		return lerp(startPos, endPos, t);
	}
};

//StructuredBuffer<float>  gRandom : register(t0);
StructuredBuffer<Sphere> gInput : register(t1);


float hit_sphere(Ray r, Sphere s)
{
	float3 oc = r.origin - s.origin(r.time);
	float a = dot(r.direction, r.direction);
	float b = dot(oc, r.direction);
	float c = dot(oc, oc) - s.radius*s.radius;
	float d = b*b - a*c;
	if (d < 0)
	{
		return -1.0;
	}
	else
	{
		return (-b - sqrt(d)) / a;
	}
}

bool my_refract(float3 v, float3 n, float ref, out float3 refracted)
{
	float dt = dot(n, v);
	float d = 1 - ref*ref*(1 - dt*dt);
	if (d > 0)
	{
		refracted = ref*(v - n*dt) - n*sqrt(d);
		return true;
	}
	return false;
}

float schlick(float cosine, float ref_idx)
{
	float r0 = (1 - ref_idx) / (1 + ref_idx);
	r0 *= r0;
	return r0 + (1 - r0)*pow(1 - cosine, 5);
}

float3 color(Ray r)
{
	float3 clr = float3(1, 1, 1);
	bool hit_anything = false;
	float t = 100000;
	float temp;
	int index;
	float3 normal;
	float ref_idx;
	float3 refracted;
	float cosine;
	float3 norm_dir;
	for (int i = 0; i < numBounces; i++)
	{
		for (int j = 0; j < numSpheres; j++)
		{
			temp = hit_sphere(r, gInput[j]);
			if (temp > 0.0001 && temp < t)
			{
				hit_anything = true;
				t = temp;
				index = j;
			}
		}
		if (!hit_anything)
		{
			float k = 0.5*(r.direction.y + 1.0);
			return clr*((1.0 - k)*white + k*light_blue);
		}
		else
		{
			hit_anything = false;
			clr *= gInput[index].material;
			normal = (r.point_at_parameter(t) - gInput[index].origin(r.time)) / gInput[index].radius;
			r.origin = r.point_at_parameter(t);
			//lambertian
			if (gInput[index].matType == 0)
			{
				r.direction = normal + random_in_unit_sphere();
			}
			else if(gInput[index].matType==1) //metal
			{
				r.direction = reflect(r.direction, normal) + gInput[index].matAdd*random_in_unit_sphere();
				if (dot(r.direction, normal) < 0)
					return float3(0.0, 0.0, 0.0);
			}
			else //dielectric
			{
				norm_dir = normalize(r.direction);
				cosine = dot(norm_dir, normal);
				if ( cosine > 0)
				{
					normal = -normal;
					ref_idx = gInput[index].matAdd;
				}
				else
				{
					ref_idx = 1.0/gInput[index].matAdd;
					cosine = -cosine;
				}
				if (all((refracted = refract(norm_dir, normal, ref_idx))==float3(0,0,0)))
				{
					r.direction = reflect(r.direction, normal);
					
				}
				else
				{
					if (random(0.0,1.0) < schlick(cosine, ref_idx))
					{
						r.direction = reflect(r.direction, normal);
					}
					else
					{
						r.direction = refracted;
					}
					
				}
			}
			t = 100000;
		}
	}

	return float3(0.0, 0.0, 0.0);
}

float3 color2(Ray r)
{
	float3 clr = float3(1.0, 1.0, 1.0);
	bool hit_anything = false;
	float t = 100000;
	float temp;
	int index;
	float3 normal;
	float ref_idx;
	float3 refracted;
	float cosine;
	float3 norm_dir;
	int r1;
	int r2;

	for (int i = 0; i < numBounces; i++)
	{
		
		for (int j = 0; j < numSpheres; j++)
		{

				temp = hit_sphere(r, gInput[j]);
				if (temp > 0.0001 && temp < t)
				{
					hit_anything = true;
					t = temp;
					index = j;
				}
			
		}
		if (!hit_anything)
		{
			float k = 0.5*(r.direction.y + 1.0);
			return clr*((1.0 - k)*white + k*light_blue);
		}
		else
		{
			hit_anything = false;
			clr *= gInput[index].material;
			r.origin = r.point_at_parameter(t);
			normal = (r.origin - gInput[index].origin(r.time)) / gInput[index].radius;
			
			//lambertian

			if (gInput[index].matType == 0)
			{
				r.direction = normal + random_in_unit_sphere();
			}
			else if (gInput[index].matType == 1) //metal
			{
				if (dot(r.direction, normal) >  0)
					return float3(1.0,0,0);
				r.direction = reflect(r.direction, normal) + gInput[index].matAdd*random_in_unit_sphere();
				
			}
			else //dielectric
			{
				norm_dir = normalize(r.direction);
				cosine = dot(norm_dir, normal);

				if (cosine > 0)
				{
					normal = -normal;
					ref_idx = gInput[index].matAdd;
				}
				else
				{
					ref_idx = 1.0 / gInput[index].matAdd;
					cosine = -cosine;
				}
				if (!any((refracted = refract(norm_dir, normal, ref_idx))))
				{
					r.direction = reflect(r.direction, normal);

				}
				else
				{
					r.direction = refracted;

				}
			}
			t = 100000;


		}
	}

	return clr/5;
}

float rand(float2 co) {
	return frac(sin(dot(co.xy, float2(12.9898, 78.233))) * 43758.5453);
}

float4 main(VS_OUTPUT input) : SV_TARGET
{
	
	//seed(randomVal*rand_1_05(float2(input.p.x,input.p.y)));
	/*float xindex = (input.p.x + 1.0) / 2.0;
	float yindex = (input.p.y + 1.0) / 2.0;
	int index = xindex*(width - 1) + width*(height - 1)*yindex;
	int seed = gRandom[index]*MAX_INT;*/
	int seed = rand(float2(input.p.x, input.p.y))*MAX_INT;
	rd.SetSeed(seed);
	Ray r;
	r.origin = origin;
	float3 clr = float3(0.0, 0.0, 0.0);
	for (int i = 0; i < numSamples; i++)
	{
		r.time = random(0.0, 1.0);
		r.direction = center+(input.p.x + random(0.0, 0.0001))*horizontal + (input.p.y + random(0.0, 0.0001))*vertical;
		clr += color2(r);
	}
	return float4(sqrt(clr / numSamples),1.0);
}