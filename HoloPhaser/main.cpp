
#include <QtCore/QCoreApplication>
#include <QFile>
#include "GSATracking.h"
#include <iostream>

int main(int argc, char *argv[])
{
	QString ampFile("../example_data/amplitude_199x199_double.raw");
	QString phaFile("../example_data/phase_199x199_double.raw");
	QString solFile("../example_data/solution_199x199_double.raw");

	QSize rawSize(199,199);
	int sz = rawSize.width() * rawSize.height();

	QFile amplitude(ampFile);
	if (!amplitude.open(QIODevice::ReadOnly))
	{
		qDebug() << "amplitude file not found" << ampFile;
		return -1;
	}

	double* amp = new double[sz];
	amplitude.read((char*)amp, sz * sizeof(amp[0]));	
	amplitude.close();


	QFile phase(phaFile);
	if (!phase.open(QIODevice::ReadOnly))
	{
		qDebug() << "phase file not found" << phaFile;
		return -1;
	}

	double* pha = new double[sz];
	phase.read((char*)pha, sz * sizeof(pha[0]));
	phase.close();

	GSATracking up;
	up.SetBranchCutMode(GSATracking::AMPLITUDETRACING); // aplitude tracking mode
	up.SetCutsDilate(7); // dilatation kernel radius
	up.SetInterpolate(true); // interpolation enabled
	up.SetWeightedInterpolation(true); // amplitude based weighted interpolation

	double* sol = up.UnWrapp(pha, amp, rawSize);


	QFile solution(solFile);
	solution.open(QIODevice::WriteOnly);
	solution.write((char*)sol, sz * sizeof(sol[0]));
	solution.close();


	delete amp;
	delete pha;

	qDebug() << "positive residues:" << up.PositiveResidueCount();
	qDebug() << "negative residues:" << up.NegativeResidueCount();
	qDebug() << "solution file created:" << solFile;

	qDebug() << endl << "press any key";
	std::cin.get();
	return 0;
}
