 // Copyright (c) 2015, Lukas Kvasnica
 // All rights reserved.
    
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are met:
    
 // * Redistributions of source code must retain the above copyright
 //   notice, this list of conditions and the following disclaimer.
 // * Redistributions in binary form must reproduce the above copyright
 //   notice, this list of conditions and the following disclaimer in the
 //   documentation and/or other materials provided with the distribution.
 // * Neither the name of the copyright holders nor the names of any
 //   contributors may be used to endorse or promote products derived
 //   from this software without specific prior written permission.
    
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 // ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 // POSSIBILITY OF SUCH DAMAGE.


#include "dilate.h"
#include <cstring>
#include <cmath>

void ImageTrasforms::Dilate(bool *src, bool *_dest, int sx, int sy, int halfSizeOfKernel)
{	
	bool inField = false;
	bool *dest = _dest;
	if(src == _dest)
	{
		inField = true;
		dest = new bool[sx * sy];
	}

	memset(dest, 0, sx * sy);
	for(int y = 0; y < sy; y++)
	{
		bool *row = src + y * sx;
		bool *newrow = dest + y * sx;
		for(int x = 0; x < sx; ++x)
		{
			if(row[x])
			{			
				for(int yid = -halfSizeOfKernel; yid < halfSizeOfKernel; ++yid)
				{
					int yi = y + yid;
					if( yi < 0 || yi >= sy)
						continue;

					bool *newnewrow = dest + yi * sx;
					int xlim = sqrt((double)halfSizeOfKernel * halfSizeOfKernel - (double)yid * yid); //kruhova maska

					for(int xi = x - xlim; xi < x + xlim; ++xi)
					{
						if(xi < 0 || xi >= sx)
							continue;						
						newnewrow[xi] = true;
					}
				}
			}
		}
	}

	if(inField)
	{
		memcpy(_dest, dest, sx * sy);
		delete[] dest;
	}
}

void ImageTrasforms::DilateFuzzyGauss(bool *src, double *dest, int sx, int sy, int halfSizeOfKernel)
{
	int kernelsize = 4 * halfSizeOfKernel * halfSizeOfKernel;
	double* kernel = new double[kernelsize];

	for(int yid = -halfSizeOfKernel; yid < halfSizeOfKernel; ++yid)
	{
		double* row = kernel + (yid + halfSizeOfKernel) * 2 * halfSizeOfKernel;
		for(int xid = -halfSizeOfKernel; xid < halfSizeOfKernel; ++xid)
		{
			double dl = 3 * sqrt((double)xid * xid + yid * yid) / halfSizeOfKernel;
			dl = - dl * dl / 2.;
			row[xid + halfSizeOfKernel] = exp(dl);
		}	
	}

	double sum = 0;
	for(int i = 0; i < kernelsize; ++i)
		sum += kernel[i];

	for(int i = 0; i < kernelsize; ++i)
		kernel[i] /= sum;

	memset(dest, 0, sx * sy * sizeof(dest[0]));

	for(int y = 0; y < sy; y++)
	{
		double *row = dest + y * sx;
		bool *srow = src + y * sx;
		for(int x = 0; x < sx; ++x)
		{
			double sum = 0;
			for(int yid = -halfSizeOfKernel; yid < halfSizeOfKernel; ++yid)
			{
				int yi = y + yid;
				if( yi < 0 || yi >= sy)
					continue;

				bool *newnewrow = src + yi * sx;
				double* krow = kernel + (yid + halfSizeOfKernel) * 2 * halfSizeOfKernel;

				for(int xid = -halfSizeOfKernel; xid < halfSizeOfKernel; ++xid)
				{
					int xi = x + xid; 
					if(xi < 0 || xi >= sx)
						continue;						

					if(newnewrow[xi])
						sum += krow[xid + halfSizeOfKernel];
				}
			}
			row[x] = 1 - sum;
		}
	}

	delete[] kernel;
}

void ImageTrasforms::Erode(bool *src, bool *_dest, int sx, int sy, int halfSizeOfKernel)
{
	bool inField = false;
	bool *dest = _dest;
	if(src == _dest)
	{
		inField = true;
		dest = new bool[sx * sy];
	}

	for(int y = 0; y < sy; y++)
	{
		bool *row = src + y * sx;
		bool *newrow = dest + y * sx;
		for(int x = 0; x < sx; ++x)
		{
			bool result = true;
			for(int yid = -halfSizeOfKernel; yid < halfSizeOfKernel; ++yid)
			{
				int yi = y + yid;
				if( yi < 0 || yi >= sy)
					continue;

				int xlim = sqrt((double)halfSizeOfKernel * halfSizeOfKernel - (double)yid * yid); //kruhova maska

				bool *newnewrow = src + yi * sx;
				for(int xi = x - xlim; xi < x + xlim; ++xi)
				{
					if(xi < 0 || xi >= sx)
						continue;

					result = result && newnewrow[xi];
				}
			}

			newrow[x] = result;
		}
	}

	if(inField)
	{
		memcpy(_dest, dest, sx * sy);
		delete[] dest;
	}
}
