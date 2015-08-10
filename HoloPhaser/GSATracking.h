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


#ifndef GOLDSTEIN_H
#define GOLDSTEIN_H

#include <QObject>
#include <QSize>
#include <QPoint>
#include <cmath>
#include <QDebug>
#include <limits>

struct point_pairs
{
	QPoint negative;
	QPoint positive;
};

//! hashovaci funkce pro QPoint a QHash
uint qHash(QPoint key); 


class GSATracking : public QObject
{
	Q_OBJECT

public:
	GSATracking(QObject* parent = 0);
	~GSATracking();

	enum MODE{ STANDARD = 0, AMPLITUDETRACING };
	enum START{ CENTRE = 0, MAXAMPLITUDE };

	MODE BranchCutMode() { return mode; }; 


	int PositiveResidueCount() { return PRcnt; };  // count of positive residues in phase data
	int NegativeResidueCount() { return NRcnt; };  // count of negative residues in phase data
	int UnwrappedPixelsCount() { return unwrappedCount; }; // count of sucessfully unwrapped pixels

	QPoint CalcStartPoint(); // position of the phase unwrapping origin 

	double* UnWrapp(double* phase, double* amplitude, QSize _size); // the function performing phase unwrapping

public slots:
	void SetBranchCutMode(MODE _mode ) { mode = _mode; }; // method of branch cut creation STANDARD - straight lines, AMPLITUDETRACING - heuristically generated by tracing anplitude
	void SetCutsDilate(int dil) { dilate = dil; }; // sets size of dilatation kernel based on size of microscope psf
	void SetInpaintIterations(int iterations) { inpaintIterations = iterations; }; // numer of inpainting iterations for missing phase filling 
	void SetWeightedInterpolation(bool interp) { weightedInterpolation =  interp; }; // enables amplitude based calculating of missing phase in unwrapped phase
	void SetInterpolate(bool _value = true) { bInterpolate = _value; }; // enables interpolation
	void SetStartPointSelection(START _value) { startpoint = _value; }; // phase unwrapping will start from CENTRE or pixel with MAXAMPLITUDE 

	void SetMaxResiduaDistance(int _distance = -1); //< maximum distance for oposite residues connection -1 = automatic

	void SetMaxReziduaDistConst(double _resconst);

signals:
	void PResidues(int _value);
	void NResidues(int _value);
	void Unwrapped(QString _value);

	void StartPoint(QString _value);
	void MaxResiduaDistanceEquation(QString _value);
	
protected:
	int maxAmplIdx;
	int maxPairLen;
	int unwrappedCount;

	double resconst;
	int dilate;
	bool weightedInterpolation;
	bool bInterpolate;
	START startpoint; 
	MODE mode;
	int inpaintIterations; 
	int NRcnt, PRcnt;

	QSize size; // image size

	double* Phase;
	double* Amplitude;
	int sz;	

	double* solution;
	bool* unwrappedPixels,* branchCuts;

	QList<QPoint> PositiveResidues;
	QList<QPoint> NegativeResidues;
	QList<point_pairs> ResiduesPairs;
	
	QList<int> TrackList;

	void UnWrap(); 
	void DetectResidua();
	void PairResidua();
	void GenerateBranchCuts();
	void InterpolateCuts();
	inline double Modulo(const double &_a, const double &_b){ return _a - floor(_a / _b) * _b; }

    void TraceAmplitude();
	void DilateCuts();

	QList<QPoint> branches; 
	QList<QPoint> trace; //< seznam bodu detekovanych minim
	void InpaintPhase();

};




#endif // GOLDSTEIN_H