//WeeklyPhysicianScheduling
#include <stdio.h>  
#include <conio.h>  
#include <stdlib.h>  
#include <time.h>  
#include <iostream>  
#include <iomanip>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <math.h> 
#include <string> 
#include <fstream> 
#include <algorithm>
using namespace std;

ILOSTLBEGIN

typedef IloArray<IloNumArray> Float2D;
typedef IloArray<Float2D> Float3D;
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
typedef IloArray<IloNumVarArray4> IloNumVarArray5;
typedef IloArray<IloNumVarArray5> IloNumVarArray6;

const bool isShowDetails = true;

const int numberOfTotalPhysicians = 16;
const int numberOfFinicPhysicians = 10;
const int numberOfDays = 7;
const int numberOfPeriods = 24;

double punishCostPerCustomerInQueue = 10.0;
double punishCostPerPhy = 2.0;
double thresholdOfCustomerInQueue = 15.0;
int minimizePhysicianInEachPeriod = 1;

const int numberOfLinearizedPieces = 10;

double slope[numberOfTotalPhysicians+1][numberOfLinearizedPieces+1] = {0};
double inter[numberOfTotalPhysicians+1][numberOfLinearizedPieces+1] = {0};
string addressOfSlopeAndInter = "separation_p10_wck.txt";

double lambda[numberOfDays][numberOfPeriods] = {
		{18.89,27.71,20.55,23.97,18.78,14.63,17.02,16.91,15.25,20.03,17.12,15.56,26.98,23.66,16.39,12.56,8.3,7.99,5.4,3.32,2.91,2.7,6.23,7.16},
		{19.07,33.57,22.99,23.65,16.34,14.78,9.07,17.08,13.89,14.85,18.04,14.48,26.81,28.14,30.04,18.04,11.12,6.59,4.61,3.15,4.72,2.15,4.27,8.25},
		{17.2,33.19,22.08,19.67,19.99,17.66,10.3,9.28,11.22,15.44,15.31,16.68,13.55,22.31,28.49,23.19,12.34,7.59,6.26,4.29,3.08,3.74,3.62,6.36},
		{23.89,41.12,11.3,20.03,17.94,11.05,14.95,19.14,20.48,23.48,27.79,28.64,25.54,26.31,21.65,11.61,6.03,4.46,2.45,2.6,1.67,2.92,1.14,3.55},
		{9.24,28.14,24.89,11.05,14.53,16.45,9.79,16.52,20.73,9.18,13.33,12.86,21.04,27.68,23.02,15.5,11.5,7.26,10.68,5.14,3.99,3.25,3.72,5.22},
		{21.55,33.27,23.31,23.25,15.87,14.53,7.25,16.05,12.29,13.39,16.49,10.96,26.24,28.06,30.32,18.32,10.86,6.44,4.43,2.63,4.68,1.62,3.72,6.3},
		{21.53,39.84,26.68,18.78,15.88,11.8,11.68,15.06,15,10.93,11.02,11.3,26.48,27.26,19.64,12.65,6.19,6.94,3.4,6.14,3.26,2.57,2.23,7.58}
};
double mu = 5.9113;
double readLambda(int d, int t) { return lambda[d - 1][t - 1]; }

void readSlopeAndInter(string filename) throw (exception){
	ifstream infile;
	try {
		infile.open(filename);
		if (!infile) {
			cout << "open file failed" << endl; system("pause");
		}

		for (int k = 1; k <= numberOfTotalPhysicians; k++) {
			for (int v = 1; v <= numberOfLinearizedPieces; v++) {
				string tempLine;
				getline(infile,tempLine);
				
				std::stringstream ss(tempLine);
				std::string temp;
				std::vector<string>elems;
				while (std::getline(ss, temp, ' ')) {
					if (!temp.empty()) {
						elems.push_back(temp);
					}
				}
				slope[k][v] = atof(elems[0].c_str());
				inter[k][v] = atof(elems[1].c_str());
			}
			string tempLine1;
			getline(infile, tempLine1);			
		}
	}
	catch (exception e) {
		e.what();
		}
}

int maxNumberOfShiftsPerDay = 2;
int lbOfshiftLength = 5;
int ubOfshiftLength = 8;
int totalMaxShiftLength = 10;
int minGapOfShifts = 2;
int startOfNigthtShift = 17;
int endOfNightShift = 24;

class Shift {
public:
	int startPeriod;
	int endPeriod;
	Shift(int a1, int a2) {
		this->startPeriod = a1;
		this->endPeriod = a2;
	}
	static bool isNotConflict(Shift s1, Shift s2) {
		if ((s1.endPeriod + minGapOfShifts < s2.startPeriod) ||
			(s2.endPeriod + minGapOfShifts < s1.startPeriod)) {
			if (((s1.endPeriod - s1.startPeriod + 1) + (s2.endPeriod - s2.startPeriod + 1)) <= totalMaxShiftLength)
				return true;
			else return false;
		}else
			return false;
	}
};

vector<Shift> shiftSets;
int initialShiftSet() {
	int countShifts = 0;
	for (int startPeriod = 1; startPeriod < startOfNigthtShift; startPeriod++) {
		for (int endPeriod = startPeriod + lbOfshiftLength - 1; endPeriod < min(startOfNigthtShift, startPeriod + ubOfshiftLength); endPeriod++) {
			Shift newShift(startPeriod,endPeriod);
			shiftSets.push_back(newShift);
		}
	}
	return shiftSets.size();
}

class DailySchedule {
	bool on_duty[numberOfPeriods+1];
public:
	void set(int period, bool state = true) { on_duty[period] = state; }
	void initial() {
		for (int i = 0; i <= numberOfPeriods; i++)on_duty[i] = false;
	}
	void setToOnDutyByShift(Shift shift) {
		for (int i = shift.startPeriod; i <= shift.endPeriod; i++)on_duty[i] = true;
	}
	bool isOnDuty(int period) { return on_duty[period]; }
};

vector<DailySchedule> scheduleSet;
int initialSchedulePlans() {
	int shiftSize = shiftSets.size();
	if (isShowDetails)cout << "the shift size(only day shift) is \t" << shiftSize << endl;

	//add day schedule with one shift
	for (int indexOfShift = 0; indexOfShift < shiftSize; indexOfShift++) {
		DailySchedule newSchedule;
		newSchedule.initial();
		newSchedule.setToOnDutyByShift(shiftSets[indexOfShift]);
		scheduleSet.push_back(newSchedule);
	}

	//add day schedule with two shifts
	for (int indexOfShift1 = 0; indexOfShift1 < shiftSize; indexOfShift1++) {
		for (int indexOfShift2 = indexOfShift1+1; indexOfShift2 < shiftSize; indexOfShift2++) {
			if (Shift::isNotConflict(shiftSets[indexOfShift1], shiftSets[indexOfShift2])) {
				DailySchedule newSchedule;
				newSchedule.initial();
				newSchedule.setToOnDutyByShift(shiftSets[indexOfShift1]);
				newSchedule.setToOnDutyByShift(shiftSets[indexOfShift2]);
				scheduleSet.push_back(newSchedule);
			}
		}
	}
	
	//add night schedule
	Shift nightShift(startOfNigthtShift, endOfNightShift);
	DailySchedule nightSchedule;
	nightSchedule.initial();
	nightSchedule.setToOnDutyByShift(nightShift);
	scheduleSet.push_back(nightSchedule);

	//add rest schedule
	DailySchedule relexSchedule;
	relexSchedule.initial();
	scheduleSet.push_back(relexSchedule);

	if (isShowDetails)cout << "the total schedules size is\t" << scheduleSet.size() << endl;
	return scheduleSet.size();
}

DailySchedule readScheduleSet(int indexOfSchedule) {
	if (indexOfSchedule > scheduleSet.size() || indexOfSchedule <= 0) { cout << "error: index\t" << indexOfSchedule << "\t over\t" << scheduleSet.size() << endl; system("pause"); }
	return scheduleSet[indexOfSchedule - 1];
}



void Cplex() {
	IloEnv env;
	IloModel model(env);

	//definition of the variables
	IloNumVarArray3 c_dmk(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++){
		c_dmk[d] = IloNumVarArray2(env, numberOfTotalPhysicians + 1);
		for (int m = 1; m <= numberOfTotalPhysicians; m++) {
			c_dmk[d][m] = IloNumVarArray(env, scheduleSet.size() + 1, 0, 1, ILOBOOL);
		}
	}

	IloNumVar T = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	IloNumVar Cost = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	IloNumVar CostOfPhy = IloNumVar(env, 0, IloInfinity, ILOFLOAT);

	IloNumVarArray2 Cost_dt(env, numberOfDays + 1); 
	for (int d = 1; d <= numberOfDays; d++) {
		Cost_dt[d] = IloNumVarArray(env, numberOfPeriods + 1, 0, INFINITY, ILOFLOAT);
	}

	IloNumVarArray2 q_dt(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		q_dt[d] = IloNumVarArray(env, numberOfPeriods + 1, 0, INFINITY, ILOFLOAT);
	}

	IloNumVarArray2 rou_dt(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		rou_dt[d] = IloNumVarArray(env, numberOfPeriods + 1, 0, 1, ILOFLOAT);
	}

	IloNumVarArray2 l_dt(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		l_dt[d] = IloNumVarArray(env, numberOfPeriods + 1, 0, INFINITY, ILOFLOAT);
	}

	IloNumVarArray2 s_dt(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		s_dt[d] = IloNumVarArray(env, numberOfPeriods + 1, 0, numberOfTotalPhysicians, ILOINT);
	}

	IloNumVarArray3 y_dtj(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		y_dtj[d] = IloNumVarArray2(env, numberOfPeriods + 1);
		for (int t = 1; t <= numberOfPeriods; t++) {
			y_dtj[d][t] = IloNumVarArray(env, numberOfTotalPhysicians + 1, 0, 1, ILOBOOL);
		}
	}

	IloNumVarArray2 onduty_dm(env, numberOfDays + 1);
	for (int d = 1; d <= numberOfDays; d++) {
		onduty_dm[d] = IloNumVarArray(env, numberOfTotalPhysicians + 1, 0, 1, ILOBOOL);
	}

	//set the objective function
	IloExpr sumOfTandC(env);
	sumOfTandC = T + Cost + CostOfPhy;
	
	//constraints (indexed according to our paper)
	//constraint 2
	IloExpr sumOfT(env);
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			sumOfT += s_dt[d][t];
		}
	}
	model.add(sumOfT == T);

	//constraint 3
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			IloExpr sum(env);
			for (int m = 1; m <= numberOfTotalPhysicians; m++) {
				for (int k = 1; k <= scheduleSet.size(); k++) {
					sum += readScheduleSet(k).isOnDuty(t) * c_dmk[d][m][k];
				}
			}
			model.add(s_dt[d][t] == sum);

		}
	}

	//constraint 4
	IloExpr sumOfCostForPhy(env);
	for (int d = 1; d <= numberOfDays; d++) {
		for (int m = numberOfFinicPhysicians + 1; m <= numberOfTotalPhysicians; m++) {
			sumOfCostForPhy += punishCostPerPhy * onduty_dm[d][m];
		}
	}
	model.add(sumOfCostForPhy == CostOfPhy);

	//constraint 5
	for (int d = 1; d <= numberOfDays; d++) {
		for (int m = numberOfFinicPhysicians + 1; m <= numberOfTotalPhysicians; m++) {
			IloExpr sum1(env);
			for (int k = 1; k <= scheduleSet.size(); k++) {
				sum1 += c_dmk[d][m][k];
			}
			model.add(sum1 == onduty_dm[d][m]);
		}
	}

	//constraint 6
	IloExpr sumOfCost(env);
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			sumOfCost += Cost_dt[d][t];
		}
	}
	model.add(sumOfCost == Cost);

	//constraint 7
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			model.add(Cost_dt[d][t] >= (q_dt[d][t] - thresholdOfCustomerInQueue) * punishCostPerCustomerInQueue);		
		}
	}

	//constraint 8
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			model.add(s_dt[d][t] >= 1);
		}
	}

	//constraint 9
	for (int d = 1; d <= numberOfDays; d++) {
		for (int m = 1; m <= numberOfFinicPhysicians; m++) {
			IloExpr sum1(env);
			for (int k = 1; k <= scheduleSet.size(); k++) {
				sum1 += c_dmk[d][m][k];
			}
			model.add(sum1 == 1);
		}
	}



	int indexOfNightShift = scheduleSet.size() - 1;
	int indexOfRelaxShift = scheduleSet.size();
	if (isShowDetails)cout << "indexOfNightShift\t" << indexOfNightShift << endl;
	if (isShowDetails)cout << "indexOfRelaxShift\t" << indexOfRelaxShift << endl;


	//constraint 10
	for (int m = 1; m <= numberOfFinicPhysicians; m++) {
		IloExpr sum2(env);
		for (int d = 1; d <= numberOfDays; d++) {
			sum2 += c_dmk[d][m][indexOfRelaxShift];
		}
		model.add(sum2 >= 0.99);
	}

	//constraint 11
	for (int d = 1; d <= numberOfDays - 1; d++) {
		for (int m = 1; m <= numberOfFinicPhysicians; m++) {
			model.add(c_dmk[d + 1][m][indexOfRelaxShift] >= c_dmk[d][m][indexOfNightShift]);
		}
	}

	//constraint 12
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 2; t <= numberOfPeriods; t++) {
			model.add(q_dt[d][t - 1] + readLambda(d,t) - l_dt[d][t] == q_dt[d][t]);
		}
	}
	for (int d = 2; d <= numberOfDays; d++) {
		model.add(q_dt[d - 1][numberOfPeriods] + readLambda(d, 1) - l_dt[d][1] == q_dt[d][1]);
	}	 
	model.add(0 + readLambda(1, 1) - l_dt[1][1] == q_dt[1][1]);
	
	//The remaining constraints correspond to constraint 12, which is linearized into several constraints 
	//constraint 14
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			for (int j = 1; j <= numberOfTotalPhysicians; j++) {
				for (int n = 1; n <= numberOfLinearizedPieces; n++) {
					model.add(q_dt[d][t] >= (
						slope[j][n] * (rou_dt[d][t] + y_dtj[d][t][j] - 1)
						+ y_dtj[d][t][j] * inter[j][n]
						));
				}
			}
		}
	}
	
	//constraint 20
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			IloExpr sum3(env);
			for (int j = 1; j <= numberOfTotalPhysicians; j++) {
				sum3 += (j * y_dtj[d][t][j]);
			}
			model.add(sum3 == s_dt[d][t]);
		}
	}

	//constraint 21
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			IloExpr sum4(env);
			for (int j = 1; j <= numberOfTotalPhysicians; j++) {
				sum4 += (y_dtj[d][t][j]);
			}
			model.add(sum4 == 1);
		}
	}

	double tempBigValueFor20A21 = mu * double(numberOfTotalPhysicians);
	//constraint 22
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			for (int j = 1; j <= numberOfTotalPhysicians; j++) {
				model.add(l_dt[d][t] <= (mu * rou_dt[d][t] * j + tempBigValueFor20A21 * (1 - y_dtj[d][t][j])));
			}
		}
	}

	//constraint 23
	for (int d = 1; d <= numberOfDays; d++) {
		for (int t = 1; t <= numberOfPeriods; t++) {
			for (int j = 1; j <= numberOfTotalPhysicians; j++) {
				model.add(l_dt[d][t] >= (mu * rou_dt[d][t] * j - tempBigValueFor20A21 * (1 - y_dtj[d][t][j])));
			}
		}
	}

	//one auxiliary constraint to eliminate symmetric solutions
	for (int d = 1; d <= numberOfDays; d++) {
		for (int m = numberOfFinicPhysicians + 1; m <= numberOfTotalPhysicians-1; m++) {
			model.add(onduty_dm[d][m] >= onduty_dm[d][m + 1]);
		}
	}

	model.add(IloMinimize(env, sumOfTandC));

	IloCplex cplex(model);
	cplex.setParam(IloCplex::TiLim, 60 * 60);
	cplex.setParam(IloCplex::Threads, 4);
	//cplex.exportModel("mdl.lp");

	cplex.solve();//Calling the solver


	if (cplex.getStatus() == IloAlgorithm::Infeasible)
	{
		env.out() << "no solution" << endl;
	}
	else
	{
		env.out() << endl << "Solution status:" << cplex.getStatus() << endl << "Solution value = " << cplex.getObjValue() << endl;
		
		
		cout << "s_dt" << endl;
		for (int d = 1; d <= numberOfDays; d++) {
			for (int t = 1; t <= numberOfPeriods; t++) {
				cout << "d:" << d << "\tt:" << t << ":\t" << cplex.getValue(s_dt[d][t]) << endl;
			}
		}
		cout << endl;
		cout << "q_dt" << endl;
		for (int d = 1; d <= numberOfDays; d++) {
			for (int t = 1; t <= numberOfPeriods; t++) {
				cout << "d:" << d << "\tt:" << t << ":\t" << cplex.getValue(q_dt[d][t]) << endl;
			}
		}
		cout << endl;

		cout << "the schedule in the first day" << endl;
		cout << "c_1-mk" << endl;
		for (int m = 1; m <= numberOfTotalPhysicians; m++) {
			for (int k = 1; k <= scheduleSet.size(); k++) {
				cout << "m:" << m << "\tk:" << k << ":\t" << cplex.getValue(c_dmk[1][m][k]) << endl;
			}
		}
		cout << endl;
	}
}


void showSlopeAndInter(int s, int n) {
	cout << s << "The slope of the " << n << "segment is : " << slope[s][n] << "\t, intercept is : " << inter[s][n] << endl;
}

void showSlopeAndInter() {
	for (int s = 1; s <= numberOfTotalPhysicians; s++) {
		cout << "The number of physician is£º" << s << endl;
		for (int v = 1; v <= numberOfLinearizedPieces; v++)
			showSlopeAndInter(s, v);
	}
}

void showSchedule(DailySchedule schedule) {
	for (int t = 1; t <= numberOfPeriods; t++)
		cout << schedule.isOnDuty(t);
	cout << endl;
}

void showSchedule() {
	for (int i = 1; i <= scheduleSet.size(); i++) {
		cout <<"indexOfSchedule:" << i << ":";
		showSchedule(scheduleSet[i - 1]);
	}
}

time_t timeS[1000] = { 0 }; time_t timeE[1000] = { 0 };
void timeStart(int i) {
	timeS[i] = clock();
}
void timeEndAndShow(int i) {
	timeE[i] = clock();
	cout << "time" << i << " is: " << (double(timeE[i]) - double(timeS[i])) / CLOCKS_PER_SEC<<" second"<<endl;
}


int main() {
	timeStart(1);//time function
	
	//initial the algorithm
	readSlopeAndInter(addressOfSlopeAndInter);
	showSlopeAndInter();
	initialShiftSet();
	initialSchedulePlans();

	//show the feasible daily schedules set W
	showSchedule();

	//call the slover
	Cplex();

	timeEndAndShow(1);//time function

	return 0;
}


