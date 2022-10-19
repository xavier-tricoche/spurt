int cint(float x)
{
	float fractpart, intpart;
	fractpart = fabs(modf (x , &intpart));

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

float CopySign(float x, float y) 
{
	if (y >= 0)
		return fabs(x);
	else
		return fabs(x) * -1;
}

bool IsBoundary(const int x, const int y, const int z, const int width, const int height, const int depth)
{
	if ((x < 1) || (y < 1) || (z < 1) || (x > width - 2) || (y > height - 2) || (z > depth - 2))
		return true;
	else return false;
}

bool IsValid(const int x, const int y, const int z, const int width, const int height, const int depth)
{
	if ((x < 0) || (y < 0) || (z < 0) || (x > width - 1) || (y > height - 1) || (z > depth - 1))
		return false;
	else return true;
}

bool InDomain(const float4 spt, const float4 minB, const float4 maxB)
{
	if ((spt.x < minB.x) || (spt.x > maxB.x)
	 || (spt.y < minB.y) || (spt.y > maxB.y)
	 || (spt.z < minB.z) || (spt.z > maxB.z))
		return false;
	return true;
}

int Coord2Addr(const int x, const int y, const int z, const int width, const int height, const int depth)
{
	return x + width * (y + height * z);
}


int3 Addr2Coord(const int i, const int width, const int height, const int depth)
{
	int x = i % width;
	int y = ((i - x) / width) % height;
	int z = (((i - x) / width) - y) / height;

	return ((int3)(x, y, z));
}

float4 Grid2Space(const float4 gpt, const int w, const int h, const int d, const float4 minB, const float4 maxB)
{
	float4 spt;

	spt.x = minB.x + (maxB.x - minB.x) * gpt.x / (w - 1.0);
	spt.y = minB.y + (maxB.y - minB.y) * gpt.y / (h - 1.0);
	spt.z = minB.z + (maxB.z - minB.z) * gpt.z / (d - 1.0);
	spt.w = 0.0;

	return spt;
}

float4 Space2GridF(const float4 spt, const int w, const int h, const int d, const float4 minB, const float4 maxB)
{
	float4 index;
	index.x = ((w - 1.0) * (spt.x - minB.x) / (maxB.x - minB.x));
	index.y = ((h - 1.0) * (spt.y - minB.y) / (maxB.y - minB.y));
	index.z = ((d - 1.0) * (spt.z - minB.z) / (maxB.z - minB.z));
	index.w = spt.w;
	return index;
}

int3 Space2GridI(const float4 spt, const int w, const int h, const int d, const float4 minB, const float4 maxB)
{
	int3 index;
	index.x = cint((w - 1.0) * (spt.x - minB.x) / (maxB.x - minB.x));
	index.y = cint((h - 1.0) * (spt.y - minB.y) / (maxB.y - minB.y));
	index.z = cint((d - 1.0) * (spt.z - minB.z) / (maxB.z - minB.z));
	return index;
}

float4 TimeInterpolate(float t0, float t1, float4 V0, float4 V1, float t)
{
	float frac = (t1 - t) / (t1 - t0);
	return V0 * frac + V1 * (1.0 - frac);
}

float4 GetFlowAt(float time, float4 spt, 
				 float time_0, __read_only image3d_t flow_0, sampler_t samp_0,
				 float time_1, __read_only image3d_t flow_1, sampler_t samp_1,
				 uint width, uint height, uint depth,
				 const float4 minB, const float4 maxB, bool* failed)
{
	if (!InDomain(spt, minB, maxB))
	{
		//*failed = true;
		//return (float3)(0.0,0.0,0.0);
	}
	
	float4 pos = (float4)(spt.x, spt.y, spt.z, 0.0);
	pos = Space2GridF(spt, width, height, depth, minB, maxB);
	float4 s0 = read_imagef(flow_0, samp_0, pos);
	float4 s1 = read_imagef(flow_1, samp_1, pos);
		
	return TimeInterpolate(time_0, time_1, s0, s1, time);
}

void Integrate(__global float4* y, float st, float et, 
			  float time_0, __read_only image3d_t flow_0, sampler_t samp_0,
			  float time_1, __read_only image3d_t flow_1, sampler_t samp_1,
			  const float4 minB, const float4 maxB,
			  uint width, uint height, uint depth)
{
	volatile float h = fabs(time_1 - time_0) / 500.0;
	float tcur = st;
	float tmp = et - tcur;

	// Set stepsize for integration in the direction from t to tout
	h = CopySign(h, tmp);

	// Step by step integration - as an endless loop over steps
	while(1) {
		// Adjust stepsize if necessary to hit the output point.
		// Look ahead two steps to avoid drastic changes in the stepsize and
		// thus lessen the impact of output points on the code.
		tmp = et - tcur;
		if (fabs(tmp) <= fabs(h))
			h = tmp;

		// Advance an approximate solution over one step of length h
		//////////////////////////////////////////////////////////////////////////
		float4 v1, v2, v3, v4;
		bool failed = false;

		v1 = h * GetFlowAt(tcur, (*y), 
						   time_0, flow_0, samp_0, time_1, flow_1, samp_1, width, height, depth, minB, maxB, &failed);
		
		
		v2 = h * GetFlowAt(tcur + 0.5f * h, (*y) + 0.5f * v1, 
						   time_0, flow_0, samp_0, time_1, flow_1, samp_1, width, height, depth, minB, maxB, &failed);

			
		float4 oy = (*y);
		if (!failed) (*y) = (*y) + v2;
		//////////////////////////////////////////////////////////////////////////

		if (fabs(tmp) <= fabs(h)) {
			// integration complete
			(*y).w = 0.0;
			
			return;
		}

		if (failed) {
			// left volume boundaries
			if (fabs(h) < 0.0000001)
			{
				(*y).w = 1.0;
				return;
			}
				
			// back track and reduce step
			*y = oy;
			h = 0.05 * h;
			tcur -= h;
		}

		tcur += h;
	}
}


__kernel void
d_advectParticles(	__global float4* flow_map, 
					uint width, uint height, uint depth,
					float4 minB4, float4 maxB4, float4 spacing4,
					float time_0, __read_only image3d_t flow_0, sampler_t samp_0,
					float time_1, __read_only image3d_t flow_1, sampler_t samp_1,
					float st, float et, uint boundary_handling, int portion
				 )

{
	uint idx = width * height * portion + get_global_id(0);
	
	if (flow_map[idx].w == 1.0)
		return;
		
	// float4 minB = (float3)(minB4.x, minB4.y, minB4.z);
	// float4 maxB = (float3)(maxB4.x, maxB4.y, maxB4.z);
	// float4 spacing = (float3)(spacing4.x, spacing4.y, spacing4.z);
		
	Integrate(&(flow_map[idx]), st, et, 
			  time_0, flow_0, samp_0,
			  time_1, flow_1, samp_1,
			  minB4, maxB4,
			  width, height, depth);
}