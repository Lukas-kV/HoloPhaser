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


#include "GSATracking.h" 
#include "mathconstants.h"
#include <QTime>
#include "dilate.h" 


uint qHash(QPoint key) { return  qHash(key.x() ^ key.y()); };

GSATracking::GSATracking(QObject* parent)
	: QObject(parent)
{	
	setObjectName("GoldsteinWithAmplitudeTracking");
	qsrand(QTime::currentTime().msec());

	unwrappedPixels = 0;
	solution = 0;
	branchCuts = 0;

	size = QSize(-1,-1);

	maxAmplIdx = 0;
	maxPairLen = -1;
	unwrappedCount = 0;


	bInterpolate =	true;
	startpoint =	CENTRE;
	resconst =		4.;

	//Amplitude tracking
	mode =					STANDARD;
	inpaintIterations =		200;
	dilate =				0;
	weightedInterpolation = false;

}

GSATracking::~GSATracking()
{

	delete[] unwrappedPixels;
	delete[] solution;
	delete[] branchCuts;
}

double* GSATracking::UnWrapp(double* phase, double* amplitude, QSize _size)
{
	if(mode == AMPLITUDETRACING && !amplitude)
	{
		qDebug() << "error amplitude missing for GOLDSTEIN unwrapping";
		return 0;
	}

	if(size != _size)
	{
		size = _size;
		sz = size.width() * size.height();

		delete[] solution;
		solution = new double[sz];
		
		delete[] branchCuts;
		branchCuts = new bool[sz];
		memset(branchCuts, 0, sz * sizeof(branchCuts[0]));

		SetMaxReziduaDistConst(resconst); 

		delete[] unwrappedPixels;
		unwrappedPixels = new bool[sz];
	}

	Phase = phase;
	Amplitude = amplitude;

	memset(unwrappedPixels, 0, sz * sizeof(unwrappedPixels[0]));
	memset(branchCuts, 0, sz * sizeof(branchCuts[0]));
	memset(solution, 0, sz * sizeof(solution[0]));

	UnWrap();

	unwrappedCount = 0;
	for(int i = 0; i < sz; ++i)
		if(unwrappedPixels[i])
			++unwrappedCount;

	emit Unwrapped(QString("%1% (%2/%3)").arg(QString::number((double)unwrappedCount / sz * 100., 'f', 2), QString::number(unwrappedCount), QString::number(sz)));

	CalcStartPoint();

	return solution;
}

void GSATracking::UnWrap()
{
	//volba pocatecniho bodu 
	//zvolim podle maxima amplitudy
	if(!Amplitude)
		return;

	DetectResidua();
    
	switch(mode)
	{
		case STANDARD:			
			    PairResidua();
				GenerateBranchCuts();			
				break;
		case AMPLITUDETRACING:
				TraceAmplitude();				
				break;
	}
	// nalezeni pocatecniho bodu navazovani faze dle polohy maxima amplitudy	
	double max = Amplitude[0];
	maxAmplIdx = 0; //sz / 2;
	for(int i = 0; i < sz; ++i)
		if(!branchCuts[i] && max < Amplitude[i])
		{
			max = Amplitude[i];
			maxAmplIdx = i;
		}

	switch(startpoint)
	{
		case 0:
			solution[sz/2] = Phase[sz/2];
			unwrappedPixels[sz/2] = true;
			TrackList.append(sz/2);
			break;

		case 1:
			solution[maxAmplIdx] = Phase[maxAmplIdx];
			unwrappedPixels[maxAmplIdx] = true;
			TrackList.append(maxAmplIdx);
			break;
	}

	int pidx = 0; // index pixelu z tracklistu
	int px = 0, py = 0; // souradnice pixelu z tracklistu
	double C = 0.; // hodnota pixelu z tracklistu

	int idx = 0; // index sousedniho pixelu k pixelu v tracklistu
	double N = 0.; // hodnota sousedniho pixelu
	double diff = 0.; // rozdil hodnot N a C;

	while(!TrackList.isEmpty())
	{
		pidx = TrackList.first();
		//pidx = TrackList.at(((double)qrand() / RAND_MAX) * TrackList.size() - 1);
		px = pidx % size.width();
		py = pidx / size.width();
		C = solution[pidx];

		// pro kazdy ze ctyr sousednich pixelu
		if(px - 1 >= 0)
		{
			idx = pidx - 1;
			if(!unwrappedPixels[idx] && !TrackList.contains(idx) && !branchCuts[idx])
			{
				// unwrapp
				N = Phase[idx];
				diff = C - N;

				if(diff > PI) 
				{
					while(C - N > PI)
						N += M_2PI;			
				}
				else
					if(diff < -PI) 
					{
						while(C - N < -PI)
							N -= M_2PI;
					}

					solution[idx] = N;
					TrackList.append(idx);
					unwrappedPixels[idx] = true;
			}
		}

		if(px + 1 < size.width())
		{
			idx = pidx + 1;
			if(!unwrappedPixels[idx] && !TrackList.contains(idx) && !branchCuts[idx])
			{
				// unwrapp
				N = Phase[idx];
				diff = C - N;

				if(diff > PI) 
				{
					while(C - N > PI)
						N += M_2PI;			
				}
				else
					if(diff < -PI) 
					{
						while(C - N < -PI)
							N -= M_2PI;
					}

					solution[idx] = N;
					TrackList.append(idx);
					unwrappedPixels[idx] = true;
			}
		}
	
		if(py - 1 >= 0)
		{
			idx = pidx - size.width();
			if(!unwrappedPixels[idx] && !TrackList.contains(idx) && !branchCuts[idx])
			{
				// unwrapp
				N = Phase[idx];
				diff = C - N;

				if(diff > PI) 
				{
					while(C - N > PI)
						N += M_2PI;			
				}
				else
					if(diff < -PI) 
					{
						while(C - N < -PI)
							N -= M_2PI;
					}

					solution[idx] = N;
					TrackList.append(idx);
					unwrappedPixels[idx] = true;
			}
		}


		if(py + 1 < size.height())
		{
			idx = pidx + size.width(); 
			if(!unwrappedPixels[idx] && !TrackList.contains(idx) && !branchCuts[idx])
			{
				// unwrapp
				N = Phase[idx];
				diff = C - N;

				if(diff > PI) 
				{
					while(C - N > PI)
						N += M_2PI;			
				}
				else
					if(diff < -PI) 
					{
						while(C - N < -PI)
							N -= M_2PI;
					}

					solution[idx] = N;
					TrackList.append(idx);
					unwrappedPixels[idx] = true;
			}
		}

		TrackList.removeFirst(); //.removeOne(pidx);
	}

	switch(mode)
	{
		case STANDARD:			
				InterpolateCuts();
				break;
		case AMPLITUDETRACING:
				DilateCuts();
				InpaintPhase();			
				break;
	}	
}

void GSATracking::DetectResidua()
{
	double active, below, belowright, right; /// hodnoty faze v okolnich pixelech
	double sum;

	int ybnd = size.height() - 1;
	int xbnd = size.width() - 1;

	double* p,* t;

	PositiveResidues.clear();
	NegativeResidues.clear();

	for(int y = 0; y < ybnd; ++y)
	{	
		p = Phase + y * size.width();

		for(int x = 0; x < xbnd; ++x)
		{
			t = p + x;
			active = *t;
			right = t[1];
			below = t[size.width()];
			belowright =t[1 + size.width()];
			
			sum = Modulo(active - below + PI, M_2PI) - PI;          //Wrap the phase differences as we loop around the 2 by 2 blocks 
			sum += Modulo(below - belowright + PI, M_2PI) - PI; 
			sum += Modulo(belowright - right + PI, M_2PI) - PI; 
			sum += Modulo(right - active + PI, M_2PI) - PI; 
			
			if(sum >= 6.) //positive residue (which should equal 2*pi)  
				PositiveResidues.append(QPoint(x,y));
			else
				if(sum <= -6) //negative residues (which should equal -2*pi) 			
					NegativeResidues.append(QPoint(x,y));
		}
	}


	NRcnt = NegativeResidues.size();
	PRcnt = PositiveResidues.size();
	
	emit PResidues(PRcnt);
	emit NResidues(NRcnt);
}

void GSATracking::PairResidua()
{
	QPoint P,N;
	int d1, d2, D;
	point_pairs pp;

	ResiduesPairs.clear();

	for(int dist = 1; dist < maxPairLen; ++dist)
	{
		
		//najdi vsechny residua mezi kterimi je vydalenost dist

		for(int a = 0; a < PositiveResidues.size(); ++a)
		{
			P = PositiveResidues.at(a);
			for(int b = 0; b < NegativeResidues.size(); ++b)
			{
				N = NegativeResidues.at(b);
				//vzdalenost definujeme jako vetsi absolutni hodnotu rozdilu souradnic dvou bodu
				d1 = abs(P.x() - N.x());
				d2 = abs(P.y() - N.y());

				if(d1 > d2)
					D = d1;
				else
					D = d2;

				if(D == dist)
				{
					PositiveResidues.removeAt(a);
					NegativeResidues.removeAt(b);
					a--;

					pp.negative = N;
					pp.positive = P;

					ResiduesPairs.append(pp);

					break;
				}			
			}		
		}	
	}

	//qDebug() << "rest negative " << NegativeResidues.size() << " positive " << PositiveResidues.size();

	if(PositiveResidues.size())
	{
		for(int a = 0; a < PositiveResidues.size(); ++a)
		{
			P = PositiveResidues.at(a);
			
			d1 = size.width() - P.x(); // vzdalenost k pravemu okraji
			d2 = size.height() - P.y(); // vydalenost ke spodnimu okraji

			bool left = (d1 > P.x());
			bool top = (d2 > P.y());			
			
			//ctyri sektory obrazo

			if(left && top)
			{
				// horni levy
				if(P.x() < P.y()) // x je blize k ose
					N = QPoint(0, P.y());
				else
					N = QPoint(P.x(), 0);
			}

			if(!left && top)
			{
				//horni pravy
				if(P.y() < d1) 
					N = QPoint(P.x(), 0);
				else
					N = QPoint(size.width() - 1, P.y());

			}

			if(left && !top)
			{
				// dolni levy
				if(P.x() < d2) 
					N = QPoint(0 , P.y());
				else
					N = QPoint(P.x(), size.height() - 1);
			}

			if(!left && !top)
			{
				//dolni pravy
				if(d1 < d2) 
					N = QPoint(size.width() - 1, P.y());
				else
					N = QPoint(P.x(), size.height() - 1);
			
			}

			pp.negative = N;
			pp.positive = P;
			ResiduesPairs.append(pp);

		}

		PositiveResidues.clear();
	}

	if(NegativeResidues.size())
	{
		for(int a = 0; a < NegativeResidues.size(); ++a)
		{
			N = NegativeResidues.at(a);
			
			d1 = size.width() - N.x(); // vzdalenost k pravemu okraji
			d2 = size.height() - N.y(); // vydalenost ke spodnimu okraji

			bool left = (d1 > N.x());
			bool top = (d2 > N.y());			
			
			//ctyri sektory obrazo

			if(left && top)
			{
				// horni levy
				if(N.x() < N.y()) // x je blize k ose
					P = QPoint(0, N.y());
				else
					P = QPoint(N.x(), 0);
			}

			if(!left && top)
			{
				//horni pravy
				if(N.y() < d1) 
					P = QPoint(N.x(), 0);
				else
					P = QPoint(size.width() - 1, N.y());

			}

			if(left && !top)
			{
				// dolni levy
				if(N.x() < d2) 
					P = QPoint(0 , N.y());
				else
					P = QPoint(N.x(), size.height() - 1);
			}

			if(!left && !top)
			{
				//dolni pravy
				if(d1 < d2) 
					P = QPoint(size.width() - 1, N.y());
				else
					P = QPoint(N.x(), size.height() - 1);
			
			}

			pp.negative = N;
			pp.positive = P;
			ResiduesPairs.append(pp);

		}

		NegativeResidues.clear();		
	}

	//qDebug() << "pairs " << ResiduesPairs.size();
}

void GSATracking::GenerateBranchCuts()
{
	QPoint P,N;
	double k,q;
	int dx, dy;

	for(int p = 0; p < ResiduesPairs.size(); ++p)
	{
		P = ResiduesPairs.at(p).positive;		
		N = ResiduesPairs.at(p).negative;

		branchCuts[P.x() + P.y() * size.width()] = true;

		if(P == N) // pokud je residum primo na hranici obrazu
			continue;

		branchCuts[N.x() + N.y() * size.width()] = true;

		dx = P.x() - N.x();
		dy = P.y() - N.y();


		if(abs(dx) > abs(dy)) 
		{
			// rovnice primky ze dvou bodu y jako funkce x kvuli vzorkovani 
			// y = kx + q

			k = (double)dy / dx;
			q = P.y() - k * P.x();
			
			if(dx > 0)
				for(int x = N.x(); x <= P.x(); ++x)
				{
					int y = k * x + q;
					branchCuts[x + y * size.width()] = true;				
				}
			else
				for(int x = P.x(); x <= N.x(); ++x)
				{
					int y = k * x + q;
					branchCuts[x + y * size.width()] = true;				
				}		
		}
		else
		{
			// rovnice primky ze dvou bodu x jako funkce y kvuli vzorkovani 
			// x = ky + q

			k = (double)dx / dy;
			q = P.x() - k * P.y();
				
			if(dy > 0)
				for(int y = N.y(); y <= P.y(); ++y)
				{
					int x = k * y + q;
					branchCuts[x + y * size.width()] = true;
				}
			else
				for(int y = P.y(); y <= N.y(); ++y)
				{
					int x = k * y + q;
					branchCuts[x + y * size.width()] = true;
				}
		}
	}
}

void GSATracking::InterpolateCuts()
{
	if(!bInterpolate)
		return;

	bool* row; 
	int x1,x2,y1,y2;
	double v1, v2, Vx, Vy;
	bool a,b;

	// interpolace (x,y)
	//
	//			(x,y1) 
	//(x1,y)	(x,y)    (x2, y)
	//			(x,y2)

	for(int y = 0; y < size.height(); ++y)
	{
		row = branchCuts + y * size.width();
	
		for(int x = 0; x < size.width(); ++x)
		{
			if(row[x]) // provest interpolaci
			{
				x1 = x2 = 0;
				y1 = y2 = 0;
				//najit nejblizsi platne body v obou osach
				
				for(int cx = x; cx < size.width(); ++cx)
					if(!row[cx])
					{
						x1 = cx;
						break;
					}

				for(int cx = x; cx >= 0; --cx)
					if(!row[cx])
					{
						x2 = cx;
						break;
					}

				for(int cy = y; cy < size.height(); ++cy)
					if(!branchCuts[x + cy * size.width()])
					{
						y1 = cy;
						break;
					}

				for(int cy = y; cy >= 0; --cy)
					if(!branchCuts[x + cy * size.width()])
					{
						y2 = cy;
						break;
					}
			
					// provereni bodu, aby interpolace fungovala					
					a = x2 != x1;
					b = y2 != y1;

					if(a && b)
					{
						//interpolace v x smeru
						v1 = solution[x1 + y * size.width()];
						v2 = solution[x2 + y * size.width()];
						Vx =  v1 + (x - x1) * (v2 - v1) / (x2 - x1);

						//interpolace v y smeru
						v1 = solution[x + y1 * size.width()];
						v2 = solution[x + y2 * size.width()];
						Vy =  v1 + (y - y1) * (v2 - v1) / (y2 - y1);

						//solution[x + y * size.width()] = (Vx + Vy) / 2; 

						if(abs(x1 - x2) < abs(y1 - y2))
							solution[x + y * size.width()] = Vx;
						else
							solution[x + y * size.width()] = Vy;
					}
					else
						if(a)
						{
							//interpolace v x smeru
							v1 = solution[x1 + y * size.width()];
							v2 = solution[x2 + y * size.width()];
							Vx =  v1 + (x - x1) * (v2 - v1) / (x2 - x1);

							solution[x + y * size.width()] = Vx;
						}
						else
							if(b)
							{
								//interpolace v y smeru
								v1 = solution[x + y1 * size.width()];
								v2 = solution[x + y2 * size.width()];
								Vy =  v1 + (y - y1) * (v2 - v1) / (y2 - y1);

								solution[x + y * size.width()] = Vy;
							}
							else
								solution[x + y * size.width()] = 0.;


			}

		}
	
	}
}

void GSATracking::TraceAmplitude()
{
	branches.clear();

    for(int i = PositiveResidues.size() - 1; i >= 0; --i )
    {
        QPoint start = PositiveResidues.at(i);

        bool run  = true; bool connected = false;
        //double m[8]; 
		QPoint p[8]; double min, stopMin;
		trace.clear();
        QPoint current = start;        
		QHash<QPoint, double> tracePoints;
        do
        {
            //if(current.x() < size.width() - 1 && current.x() > 0 && current.y() < size.height() - 1 && current.y() > 0)
            {
                trace.append(current);
                //branchCuts[current.y() * size.width() + current.x()] = true;
                //double c = Amplitude[current.y() * size.width() + current.x()];
                p[0] = QPoint(current.x(), current.y() - 1);
                p[1] = QPoint(current.x() + 1, current.y() - 1);
                p[2] = QPoint(current.x() + 1, current.y());
                p[3] = QPoint(current.x() + 1, current.y() + 1);
                p[4] = QPoint(current.x(), current.y() + 1);
                p[5] = QPoint(current.x() - 1, current.y() - 1);
                p[6] = QPoint(current.x() - 1, current.y() + 1);
                p[7] = QPoint(current.x() - 1, current.y());

                //m[0] = Amplitude[p[0].y() * size.width() + p[0].x()];
                //m[1] = Amplitude[p[1].y() * size.width() + p[1].x()];
                //m[2] = Amplitude[p[2].y() * size.width() + p[2].x()];
                //m[3] = Amplitude[p[3].y() * size.width() + p[3].x()];
                //m[4] = Amplitude[p[4].y() * size.width() + p[4].x()];
                //m[5] = Amplitude[p[5].y() * size.width() + p[5].x()];
                //m[6] = Amplitude[p[6].y() * size.width() + p[6].x()];
                //m[7] = Amplitude[p[7].y() * size.width() + p[7].x()];


				//! presunuti vsech potencialnich minim do tracepoints
                for(int j = 0; j < 8; ++j)
                {
					if(!trace.contains(p[j]) && p[j].x() > -1 && p[j].y() > -1 && p[j].y() < size.height() && p[j].x() < size.width())
                    {
						tracePoints[p[j]] = Amplitude[p[j].y() * size.width() + p[j].x()];
                    }
                }
				
				// vyhledani minimalni hodnoty

				QHash<QPoint, double>::const_iterator it = tracePoints.constBegin();
				min = it.value(); 
				current = it.key();

				while (it != tracePoints.constEnd()) 
				{
					if(it.value() < min)
					{
						min = it.value();
						current = it.key();
					}
					++it;
				}

				tracePoints.remove(current);

				if(!connected)
				{
					for(int j = 0; j < 8; ++j)
					{
						if(NegativeResidues.contains(p[j]))
						{
							NegativeResidues.removeOne(p[j]);
							PositiveResidues.removeAt(i);
							trace.append(p[j]);
							//qDebug() << "propojeno";
							connected = true;
							stopMin = Amplitude[p[j].y() * size.width() + p[j].x()];
							//run = false;
							break;
						}
					}
				}
				else
				{
					if(min > stopMin)
						run = false;
				}
            }

        } while(run && trace.size() < 1000);

		branches.append(trace);

        for(int j = 0; j < trace.size(); ++j)
            branchCuts[trace.at(j).y() * size.width() + trace.at(j).x()] = true;

		if(NegativeResidues.size() && PositiveResidues.size())
        {
            QPoint r = NegativeResidues.at(0);
            NegativeResidues[0] = PositiveResidues.at(0);
            PositiveResidues[0] = r;
        }
    }
}

void GSATracking::DilateCuts()
{
	if(dilate)
	{
		ImageTrasforms::Dilate(branchCuts, branchCuts, size.width(), size.height(), dilate);
	}
}

void GSATracking::InpaintPhase()
{
	if(!bInterpolate)
		return;

	QList<QPoint> orig = branches;
	branches.clear();
	for(int y = 0; y < size.height(); ++y)
	{
		bool* row = branchCuts + y * size.width();
		for(int x = 0; x < size.width(); ++x)
		{
			if(row[x])
				branches.append(QPoint(x,y));
		}
	}
	
	double* original = new double[branches.size()];
	double* weights = new double[branches.size()];
	int as = 0;
	foreach(QPoint p, branches)
	{
		int idx = p.y() * size.width() + p.x(); 
		original[as] = solution[idx];
		solution[idx] = 0;

		if(orig.contains(p))
			weights[as] = 0;
		else
			weights[as] = Amplitude[idx];
		++as;
	}

	double max = weights[0], min = weights[0];	
	for(int i = 0; i < branches.size(); ++i)
	{
		if(max < weights[i])
			max = weights[i];
		else
			if(min > weights[i])
				min = weights[i];
	}
	//qDebug() << max << min;

	if(max != min)
	for(int i = 0; i < branches.size(); ++i)
	{
		weights[i] =  (weights[i] - min) / (max - min);
	}

	QPoint current, p[8];
	double ph[8];
	
	double* values = new double[branches.size()];

	for(int iterate = 0; iterate < inpaintIterations; ++iterate)
	{
		for(int j = 0; j < branches.size(); ++j)
		{
			current = branches.at(j);

			if(current.x() < size.width() - 1 && current.x() > 0 && current.y() < size.height() - 1 && current.y() > 0)
			{
				p[0] = QPoint(current.x(), current.y() - 1);
				p[1] = QPoint(current.x() + 1, current.y() - 1);
				p[2] = QPoint(current.x() + 1, current.y());
				p[3] = QPoint(current.x() + 1, current.y() + 1);
				p[4] = QPoint(current.x(), current.y() + 1);
				p[5] = QPoint(current.x() - 1, current.y() - 1);
				p[6] = QPoint(current.x() - 1, current.y() + 1);
				p[7] = QPoint(current.x() - 1, current.y());


				ph[0] = solution[p[0].y() * size.width() + p[0].x()];
				ph[1] = solution[p[1].y() * size.width() + p[1].x()];
				ph[2] = solution[p[2].y() * size.width() + p[2].x()];
				ph[3] = solution[p[3].y() * size.width() + p[3].x()];
				ph[4] = solution[p[4].y() * size.width() + p[4].x()];
				ph[5] = solution[p[5].y() * size.width() + p[5].x()];
				ph[6] = solution[p[6].y() * size.width() + p[6].x()];
				ph[7] = solution[p[7].y() * size.width() + p[7].x()];

				values[j] = 0.125 * ( ph[0] + ph[1] + ph[2] + ph[3] + ph[4] + ph[5] + ph[6] + ph[7] );	
				//values[j] = 0.176765 * ( ph[0] + ph[2] + ph[4] + ph[7] )  +  0.073235 * ( ph[1] +  ph[3] +  ph[5] + ph[6] );	
			}
			else
				values[j] = 0; //TODO nahradit nejakou rozumnou hodnotou
		}

		if(weightedInterpolation)
			for(int j = 0; j < branches.size(); ++j)
			{							
				solution[branches.at(j).y() * size.width() + branches.at(j).x()] = (values[j] * (1 - weights[j]) + original[j] * weights[j]); 
			}
		else
			for(int j = 0; j < branches.size(); ++j)
			{							
				solution[branches.at(j).y() * size.width() + branches.at(j).x()] = values[j]; 
			}

	}

	delete[] original;
	delete[] values;
	delete[] weights;
}

void GSATracking::SetMaxResiduaDistance(int _distance)
{
	if(_distance == -1)
	{
		if(size.width() < size.height())
			maxPairLen = size.width() / 4;      /// maximalni detekcni polomer odvozeny od velikosti obrazu faze
		else
			maxPairLen = size.height() / 4;
	}
	else
		maxPairLen = _distance;
}

void GSATracking::SetMaxReziduaDistConst(double _resconst) 
{
	resconst = _resconst; 

	int max_dist;
	if(size.width() < size.height())
		max_dist = size.width() / resconst;     /// maximalni detekcni polomer odvozeny od velikosti obrazu faze
	else
		max_dist = size.height() / resconst;

	SetMaxResiduaDistance(max_dist);

	emit MaxResiduaDistanceEquation(QString("min(%1,%2) / %3 = %4 px").arg(QString::number(size.width()), QString::number(size.height()),QString::number(resconst), QString::number(max_dist)));
}

QPoint GSATracking::CalcStartPoint() 
{ 	
	QPoint _value = QPoint(maxAmplIdx % size.width(), maxAmplIdx / size.width()); 

	emit StartPoint(QString("(%1, %2)").arg(QString::number(_value.x()), QString::number(_value.y())));

	return _value;
}
