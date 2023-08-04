package com.ben;

public class params {
	public static final double EPS = 1e-2;
	public static int N = 12;  // the number of total physicians (including the seconded physicians)
	public static int IndexOfFinicPhysician = 10; // the number of fever clinic physicians
	public static int G = 2;  //  the minimum interval between two shifts
	public static int D = 7;  //  the number of days
	public static int T = 24;  //  the number of time periods
	public static int Lr = 15; //  the maximum queue length threshold
	public static int S_LB = 5; //  the minimum shift length
	public static int S_UB = 8; //  the maximum shift length
	public static int S_MaxForOnePhy = 10;// the maximum working time of a physician per day
	public static int DELTA = 1; //  the length of one time period
	public static int M = 8000; //  the big M
	public static double MU = 5.9113; //  the service rate
	public static int PARTITIONS = 10; //  the number of partitions in piecewise linearization
	public static int TRUST_SCOPE = 4;   // the trust region scope
	public static int LEN_NIGHT_SCH = 8; //  the length of night shift
	public static String PWL_PARAMS_PATH = "data/separation_p10.txt";   //  the path of piecewise linearization parameters
	public static String PWL_PARAMS_PATH_wck = "data/separation_p10_wck.txt";   //  the path of piecewise linearization parameters
	public static double PunishCost = 10;// the cost of violating the queue length constraint
	public static double PuhishCostForAddedPhysician = 2; // the cost of seconded a physician
	public static int MaxIteration = 2000;
	public static boolean isHardContrOfQueueLength = false;
	public static int N_2 = 8;
	public static int minOfStaffig = 1;
	public static double[][] LAM = {
			{18.89,27.71,20.55,23.97,18.78,14.63,17.02,16.91,15.25,20.03,17.12,15.56,26.98,23.66,16.39,12.56,8.3,7.99,5.4,3.32,2.91,2.7,6.23,7.16},
			{19.07,33.57,22.99,23.65,16.34,14.78,9.07,17.08,13.89,14.85,18.04,14.48,26.81,28.14,30.04,18.04,11.12,6.59,4.61,3.15,4.72,2.15,4.27,8.25},
			{17.2,33.19,22.08,19.67,19.99,17.66,10.3,9.28,11.22,15.44,15.31,16.68,13.55,22.31,28.49,23.19,12.34,7.59,6.26,4.29,3.08,3.74,3.62,6.36},
			{23.89,41.12,11.3,20.03,17.94,11.05,14.95,19.14,20.48,23.48,27.79,28.64,25.54,26.31,21.65,11.61,6.03,4.46,2.45,2.6,1.67,2.92,1.14,3.55},
			{9.24,28.14,24.89,11.05,14.53,16.45,9.79,16.52,20.73,9.18,13.33,12.86,21.04,27.68,23.02,15.5,11.5,7.26,10.68,5.14,3.99,3.25,3.72,5.22},
			{21.55,33.27,23.31,23.25,15.87,14.53,7.25,16.05,12.29,13.39,16.49,10.96,26.24,28.06,30.32,18.32,10.86,6.44,4.43,2.63,4.68,1.62,3.72,6.3},
			{21.53,39.84,26.68,18.78,15.88,11.8,11.68,15.06,15,10.93,11.02,11.3,26.48,27.26,19.64,12.65,6.19,6.94,3.4,6.14,3.26,2.57,2.23,7.58},
	};
}
