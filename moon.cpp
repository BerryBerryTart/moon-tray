// http://www.jgiesen.de/elevazmoon/basics/index.htm
// https://astronomy.stackexchange.com/questions/51505/calculate-moon-illumination-given-moon-age

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>

const double RAD = M_PI / 180.0;
const double DEG = 180 / M_PI;
const int DEBUG = 0;
const int DEBUG_RISE_SET = 0;
const double LATITUDE = 53.400002;
const double LONGITUDE = 2.983333;

struct Moon
{
  double perc;
  std::string phase;
  double azimuth;
  double altitude;
  std::string transit;
  std::string rise;
  std::string set;
};

struct MoonPos
{
  double azimuth;
  double altitude;
  double i;
  double parallax;
  double RA;
  double Decl;
  double RA_p1;   // minus a day
  double Decl_p1; // minus a day
  double RA_m1;   // plus a day
  double Decl_m1; // plus a day
};

struct MoonRiseSet
{
  std::string transit;
  std::string rise;
  std::string set;
};

struct SumLR
{
  double L;
  double R;
};

struct JDTime
{
  double JDNow;
  double JDT;
  double JDNow_Midnight;
  double JDT_Midnight;
  double Sidereal_Midnight;
};

double constrain(double d)
{
  double t = std::fmod(d, 360.0);
  if (t < 0)
  {
    t += 360;
  }
  return t;
}

double constrain180(double d)
{
  double t = std::fmod(d, 180.0);
  if (t < 0)
  {
    t += 180;
  }
  return t;
}

double constrain1(double d)
{
  double t = d;
  if (t < 0)
  {
    t += 1;
  }
  return t;
}

// time in Julian Centuries from Epick J2000.0
JDTime julianDate(int offset = 0)
{
  JDTime jd;
  std::time_t unixTime = std::time(nullptr);
  double JDNow = ((double)unixTime / 86400.0) + 2440587.5;

  // jgiesen.de
  // JDNow = 2448396.04167;

  // Astronomical Algo
  // JDNow = 2448724.5;
  const double JD2000 = 2451545.0;
  const double T = (JDNow - JD2000) / 36525.0;

  // get sidereal time at 0h

  // get local day, month and year
  struct std::tm *t = localtime(&unixTime);

  int day = t->tm_mday;
  int month = t->tm_mon + 1;
  int year = t->tm_year + 1900;

  std::string day_str = std::to_string(t->tm_mday);
  if (day < 10)
  {
    day_str = "0" + day_str;
  }
  std::string month_str = std::to_string(t->tm_mon + 1);
  if (month < 10)
  {
    month_str = "0" + month_str;
  }
  std::string year_str = std::to_string(t->tm_year + 1900);

  //  Initialize a string stream
  std::string dateString = year_str + "-" + month_str + "-" + day_str + " 00:00:00";
  std::istringstream ss(dateString);

  std::tm midnightToday = {};
  ss >> std::get_time(&midnightToday, "%Y-%m-%d %H:%M:%S");

  if (ss.fail())
  {
    printf("UH OH!!! PARSE TIME FAILED????????\n");
  }

  double unixMidnight = (double)std::mktime(&midnightToday) + offset * 86400.0;
  // printf("mn: %lf\n", unixMidnight);
  double JDMidnight = (unixMidnight / 86400.0) + 2440587.5;
  // printf("JDM: %lf\n", JDMidnight);
  const double T_0 = (JDMidnight - JD2000) / 36525.0;
  double THETA_0 = 100.46061837 + 36000.770053608 * T_0 + 0.000387933 * std::pow(T_0, 2) - (1 / 38710000) * std::pow(T_0, 3);

  jd.JDNow = JDNow;
  jd.JDT = T;
  jd.JDNow_Midnight = JDMidnight;
  jd.JDT_Midnight = T_0;
  jd.Sidereal_Midnight = constrain(THETA_0);
  return jd;
}

std::string timeDecToString(double d)
{
  double rawTime = d * 24;
  int hour = (int)rawTime;
  double whole;
  double minute_dec = std::modf(rawTime, &whole);
  int minute = minute_dec * 60;
  std::string minute_str = std::to_string(minute);
  if (minute > 0 && minute < 10)
  {
    minute_str = "0" + minute_str;
  }
  return std::to_string(hour) + ":" + minute_str;
}

double moonPercent(double i)
{
  double k = (1 + std::cos(i)) / 2;
  return k;
}

std::string getPhase(double i, double perc)
{
  double k = (std::sin(i)) / 2.0;
  std::string phase = k > 0 ? "Waxing" : "Waning";
  if (perc > 0 && perc < 0.05)
  {
    return "New Moon";
  }
  else if (perc >= 0.05 && perc < 0.45)
  {
    return phase += " Crescent";
  }
  else if (perc >= 0.45 && perc < 0.55)
  {
    if (k > 0)
    {
      return "First Quarter";
    }
    else
    {
      return "Last Quarter";
    }
  }
  else if (perc >= 0.55 && perc < 0.95)
  {
    return phase += " Gibbous";
  }
  else if (perc >= 0.95)
  {
    return "Full Moon";
  }
  return phase;
}

SumLR getLR(double D, double M, double Mp, double F, double E)
{
  // multiples [4], sine coef, cosine coef
  std::vector<std::vector<double>> LRTable = {
      {0, 0, 1, 0, 6288774, -20905355},
      {2, 0, -1, 0, 1274027, -3699111},
      {2, 0, 0, 0, 658314, -2955968},
      {0, 0, 2, 0, 213618, -569925},
      {0, 1, 0, 0, -185116, 48888},

      {0, 0, 0, 2, -114332, -3149},
      {2, 0, -2, 0, 58793, 246158},
      {2, -1, -1, 0, 57066, -152138},
      {2, 0, 1, 0, 53322, -170733},
      {2, -1, 0, 0, 45758, -204586},

      {0, 1, -1, 0, -40923, -129620},
      {1, 0, 0, 0, -34720, 108743},
      {0, 1, 1, 0, -30383, 104755},
      {2, 0, 0, -2, 15327, 10321},
      {0, 0, 1, 2, -12528, 0},

      {0, 0, 1, -2, 10980, 79661},
      {4, 0, -1, 0, 10675, -34782},
      {0, 0, 3, 0, 10034, -23210},
      {4, 0, -2, 0, 8548, -21636},
      {2, 1, -1, 0, -7888, 24208},

      {2, 1, 0, 0, -6766, 30824},
      {1, 0, -1, 0, -5163, -8379},
      {1, 1, 0, 0, 4987, -16675},
      {2, -1, 1, 0, 4036, -12831},
      {2, 0, 2, 0, 3994, -10445},

      {4, 0, 0, 0, 3861, -11650},
      {2, 0, -3, 0, 3665, 14403},
      {0, 1, -2, 0, -2689, -7003},
      {2, 0, -1, 2, -2602, 0},
      {2, -1, -2, 0, 2390, 10056},

      {1, 0, 1, 0, -2348, 6322},
      {2, -2, 0, 0, 2236, -9884},

      // next page
      {0, 1, 2, 0, -2120, 5751},
      {0, 2, 0, 0, -2069, 0},
      {2, -2, -1, 0, 2048, -4950},
      {2, 0, 1, -2, -1773, 4130},
      {2, 0, 0, 2, -1595, 0},

      {4, -1, -1, 0, 1215, -3958},
      {0, 0, 2, 2, -1110, 0},
      {3, 0, -1, 0, -892, 3258},
      {2, 1, 1, 0, -810, 2616},
      {4, -1, -2, 0, 759, -1897},

      {0, 2, -1, 0, -713, -2117},
      {2, 2, -1, 0, -700, 2354},
      {2, 1, -2, 0, 691, 0},
      {2, -1, 0, -2, 596, 0},
      {4, 0, 1, 0, 549, -1423},

      {0, 0, 4, 0, 537, -1117},
      {4, -1, 0, 0, 520, -1571},
      {1, 0, -2, 0, -487, -1739},
      {2, 1, 0, -2, -399, 0},
      {0, 0, 2, -2, -381, -4421},

      {1, 1, 1, 0, 351, 0},
      {3, 0, -2, 0, -340, 0},
      {4, 0, -3, 0, 330, 0},
      {2, -1, 2, 0, 327, 0},
      {0, 2, 1, 0, -323, 1165},

      {1, 1, -1, 0, 299, 0},
      {2, 0, 3, 0, 294, 0},
      {2, 0, -1, -2, 0, 8752}};

  int i = 0;
  double L = 0;
  double R = 0;
  for (i = 0; i < LRTable.size(); i++)
  {
    std::vector<double> row = LRTable.at(i);
    const double Mrow = row.at(1);
    double Lcoef = 0.0;
    double Rcoef = 0.0;
    if (Mrow == 2 || Mrow == -2)
    {
      Lcoef = row.at(4) * std::pow(E, 2);
      Rcoef = row.at(5) * std::pow(E, 2);
    }
    else if (Mrow == 1 || Mrow == -1)
    {
      Lcoef = row.at(4) * E;
      Rcoef = row.at(5) * E;
    }
    else
    {
      Lcoef = row.at(4);
      Rcoef = row.at(5);
    }
    L += Lcoef * std::sin(row.at(0) * D + row.at(1) * M + row.at(2) * Mp + row.at(3) * F);
    R += Rcoef * std::cos(row.at(0) * D + row.at(1) * M + row.at(2) * Mp + row.at(3) * F);
  }

  SumLR sum;
  sum.L = L;
  sum.R = R;
  return sum;
}

double getB(double D, double M, double Mp, double F, double E)
{
  std::vector<std::vector<double>> BTable = {
      {0, 0, 0, 1, 5128122},
      {0, 0, 1, 1, 280602},
      {0, 0, 1, -1, 277693},
      {2, 0, 0, -1, 173237},
      {2, 0, -1, 1, 55413},

      {2, 0, -1, -1, 46271},
      {2, 0, 0, 1, 32573},
      {0, 0, 2, 1, 17198},
      {2, 0, 1, -1, 9266},
      {0, 0, 2, -1, 8822},

      {2, -1, 0, -1, 8216},
      {2, 0, -2, -1, 4324},
      {2, 0, 1, 1, 4200},
      {2, 1, 0, -1, -3359},
      {2, -1, -1, 1, 2463},

      {2, -1, 0, 1, 2211},
      {2, -1, -1, -1, 2065},
      {0, 1, -1, -1, -1870},
      {4, 0, -1, -1, 1828},
      {0, 1, 0, 1, -1794},

      {0, 0, 0, 3, -1749},
      {0, 1, -1, 1, -1565},
      {1, 0, 0, 1, -1491},
      {0, 1, 1, 1, -1475},
      {0, 1, 1, -1, -1410},

      {0, 1, 0, -1, -1344},
      {1, 0, 0, -1, -1335},
      {0, 0, 3, 1, 1107},
      {4, 0, 0, -1, 1021},
      {4, 0, -1, 1, 833},

      // next page
      {0, 0, 1, -3, 777},
      {4, 0, -2, 1, 671},
      {2, 0, 0, -3, 607},
      {2, 0, 2, -1, 596},
      {2, -1, 1, -1, 491},

      {2, 0, -2, 1, -451},
      {0, 0, 3, -1, 439},
      {2, 0, 2, 1, 422},
      {2, 0, -3, -1, 421},
      {2, 1, -1, 1, -366},

      {2, 1, 0, 1, -351},
      {4, 0, 0, 1, 331},
      {2, -1, 1, 1, 315},
      {2, -2, 0, -1, 302},
      {0, 0, 1, 3, -283},

      {2, 1, 1, -1, -229},
      {1, 1, 0, -1, 223},
      {1, 1, 0, 1, 223},
      {0, 1, -2, -1, -220},
      {2, 1, -1, -1, -220},

      {1, 0, 1, 1, -185},
      {2, -1, -2, -1, 181},
      {0, 1, 2, 1, -177},
      {4, 0, -2, -1, 176},
      {4, -1, -1, -1, 166},

      {1, 0, 1, -1, -164},
      {4, 0, 1, -1, 132},
      {1, 0, -1, -1, -119},
      {4, -1, 0, -1, 115},
      {2, -2, 0, 1, 107}};

  int i = 0;
  double B = 0;
  for (i = 0; i < BTable.size(); i++)
  {
    std::vector<double> row = BTable.at(i);
    const double Mrow = row.at(1);
    double Bcoef = 0.0;
    if (Mrow == 2 || Mrow == -2)
    {
      Bcoef = row.at(4) * std::pow(E, 2);
    }
    else if (Mrow == 1 || Mrow == -1)
    {
      Bcoef = row.at(4) * E;
    }
    else
    {
      Bcoef = row.at(4);
    }
    B += Bcoef * std::sin(row.at(0) * D + row.at(1) * M + row.at(2) * Mp + row.at(3) * F);
  }
  return B;
}

MoonPos getMoonPos(JDTime jd_time, int useMidnight = 0)
{
  MoonPos pos;
  const JDTime JD = jd_time;
  double T = JD.JDT;
  double JDNow = JD.JDNow;

  if (useMidnight)
  {
    T = JD.JDT_Midnight;
    JDNow = JD.JDNow_Midnight;
  }

  double Lp = 218.3164477 + 481267.881234 * T - 0.0015786 * std::pow(T, 2) + (1 / 538841.0) * std::pow(T, 3) - (1 / 65194000.0) * std::pow(T, 4);
  double D = 297.8501921 + 445267.1114034 * T - 0.0018819 * std::pow(T, 2) + (1.0 / 545868.0) * std::pow(T, 3) - (1 / 113065000.0) * std::pow(T, 4);
  double M = 357.5291092 + 35999.0502909 * T - 0.0001536 * std::pow(T, 2) + (1.0 / 24490000.0) * std::pow(T, 3);
  double Mp = 134.9633964 + 477198.8675055 * T + 0.0087414 * std::pow(T, 2) + (1.0 / 69699.0) * std::pow(T, 3) - (1 / 14712000.0) * std::pow(T, 4);
  double F = 93.2720950 + 483202.0175233 * T - 0.0036539 * std::pow(T, 2) - (1 / 3526000.0) * std::pow(T, 3) + (1 / 863310000.0) * std::pow(T, 4);
  double E = 1 - 0.002516 * T - 0.0000074 * std::pow(T, 2);

  double A1 = 119.75 + 131.849 * T;
  double A2 = 53.09 + 479264.290 * T;
  double A3 = 313.45 + 481266.484 * T;

  // clamp values to make things easier
  Lp = constrain(Lp);
  D = constrain(D);
  M = constrain(M);
  Mp = constrain(Mp);
  F = constrain(F);
  E = constrain(E);
  A1 = constrain(A1);
  A2 = constrain(A2);
  A3 = constrain(A3);

  // used for moon phase
  const double Di = D * RAD;
  const double Mpi = Mp * RAD;
  const double Mi = M * RAD;
  double i = 180.0 - Di * DEG - 6.289 * std::sin(Mpi) + 2.1 * std::sin(Mi) - 1.274 * std::sin(2 * Di - Mpi) - 0.658 * std::sin(2 * Di) - 0.214 * std::sin(2 * Mpi) - 0.11 * std::sin(Di);
  i = constrain(i) * RAD;

  SumLR sumLR = getLR(D * RAD, M * RAD, Mp * RAD, F * RAD, E);
  sumLR.L += 3958 * std::sin(A1 * RAD) + 1962 * std::sin(Lp * RAD - F * RAD) + 318 * std::sin(A2 * RAD);

  double sumB = getB(D * RAD, M * RAD, Mp * RAD, F * RAD, E);
  sumB += -2235 * std::sin(Lp * RAD) + 382 * std::sin(A3 * RAD) + 175 * std::sin(A1 * RAD - F * RAD) + 175 * std::sin(A1 * RAD + F * RAD) + 127 * std::sin(Lp * RAD - Mp * RAD) - 115 * std::sin(Lp * RAD + Mp * RAD);

  // Moon Distance
  double DELTA = 385000 + sumLR.R / 1000;

  // compute the ecliptic latitude B and the longitude L
  double LAMBDA = Lp + (sumLR.L / 1000000);
  double BETA = (sumB / 1000000);

  // obliquity of ecliptic
  double EPSILON = 23.439279 - 0.01301 * T - (0.0001831 / 3600.0) * std::pow(T, 2) + (0.00200340 / 3600) * std::pow(T, 3);

  // right ascension
  double RA = std::atan2(std::sin(LAMBDA * RAD) * std::cos(EPSILON * RAD) - std::tan(BETA * RAD) * std::sin(EPSILON * RAD), std::cos(LAMBDA * RAD)) * DEG;
  // declination delta
  double Decl = std::asin(std::sin(BETA * RAD) * std::cos(EPSILON * RAD) + std::cos(BETA * RAD) * std::sin(EPSILON * RAD) * std::sin(LAMBDA * RAD)) * DEG;

  // if atan is negative, rectify
  // if (RA < 0)
  // {
  //   RA += 180;
  // }

  RA = constrain(RA);

  // sidereal time at greenwich
  double THETA_0_G = 280.46061837 + 360.98564736629 * (JDNow - 2451545.0) + 0.000387933 * std::pow(T, 2) - (1 / 38710000) * std::pow(T, 3);

  // local sideral time at longitude
  double THETA = THETA_0_G - LONGITUDE;
  // double THETA = THETA_0 + 10;

  // local hour angle
  double TAU = THETA - RA;

  THETA = constrain(THETA);
  TAU = constrain(TAU);

  // my latitude;
  double PHI = LATITUDE;
  // double PHI = 50;

  double Alt = std::asin(std::sin(PHI * RAD) * std::sin(Decl * RAD) + std::cos(PHI * RAD) * std::cos(Decl * RAD) * std::cos(TAU * RAD)) * DEG;
  double Az = std::atan2(std::sin(TAU * RAD), std::cos(TAU * RAD) * std::sin(PHI * RAD) - std::tan(Decl * RAD) * std::cos(PHI * RAD)) * DEG + 180;

  // if (Az < 0)
  // {
  //   Az += 360;
  // }
  Az = constrain(Az);

  // parallax P
  double horP = 8.794 / (DELTA / 149597870);
  double P = std::asin(std::cos(Alt * RAD) * std::sin((horP / 3600) * RAD)) * DEG;

  // bro please work
  if (DEBUG)
  {
    printf("T: %lf, Lp: %lf, D: %lf, M: %lf, Mp: %lf, F: %lf, E: %lf\n", T, Lp, D, M, Mp, F, E);
    printf("A1: %lf, A2: %lf, A3: %lf\n", A1, A2, A3);
    printf("Sigma L: %lf, Sigma R: %lf, Sigma B: %lf\n", sumLR.L, sumLR.R, sumB);
    printf("Lambda: %lf, Beta: %lf, Epsilon: %lf\n", LAMBDA, BETA, EPSILON);
    printf("RA: %lf, Decl: %lf\n", RA, Decl);
    printf("Theta: %lf, Tau: %lf\n", THETA, TAU);
    printf("Distance: %lf, Parallax: %lf\n", DELTA, P);
  }

  pos.altitude = Alt - P;
  pos.azimuth = Az;
  pos.i = i;
  pos.RA = RA;
  pos.Decl = Decl;
  pos.parallax = P;
  return pos;
}

double getDeltaT()
{
  std::time_t unixTime = std::time(nullptr);
  // get local day, month and year
  struct std::tm *t = localtime(&unixTime);

  double year = (double)t->tm_year + 1900.0;

  double dt = (year - 2000.0) / 100.0;
  return 102 + 102 * dt + 25.3 * std::pow(dt, 2) + 0.37 * (year - 2100);
}

double interpolate(double y1, double y2, double y3, double n)
{
  double a = y2 - y1;
  double b = y3 - y2;
  double c = b - a;
  return y2 + (n / 2.0) * (a + b + n * c);
}

MoonRiseSet getMoonRiseSet(MoonPos pos, JDTime jd_time)
{
  MoonRiseSet result;
  // le constants have arrived
  const double THETA_0 = jd_time.Sidereal_Midnight;
  const double P = pos.parallax;
  const double RA_1 = pos.RA_m1;
  const double Decl_1 = pos.Decl_m1;
  const double RA_2 = pos.RA;
  const double Decl_2 = pos.Decl;
  const double RA_3 = pos.RA_p1;
  const double Decl_3 = pos.Decl_p1;

  // geometric altitude
  double h_0 = 0.7275 * P - (34.0 / 60.0);

  double H_0 = std::acos((std::sin(h_0 * RAD) - std::sin(LATITUDE * RAD) * std::sin(Decl_2 * RAD)) / (std::cos(LATITUDE * RAD) * std::cos(Decl_2 * RAD))) * DEG;

  H_0 = constrain180(H_0);

  // transit
  double m0 = (RA_2 + LONGITUDE - THETA_0) / 360.0;
  m0 = constrain1(m0);

  // rise
  double m1 = m0 - (H_0 / 360.0);
  m1 = constrain1(m1);

  // set
  double m2 = m0 + (H_0 / 360.0);
  m2 = constrain1(m2);

  // corrections
  double delta_t = getDeltaT();

  double theta_0_m0 = THETA_0 + 360.985657 * m0;
  double theta_0_m1 = THETA_0 + 360.985657 * m1;
  double theta_0_m2 = THETA_0 + 360.985657 * m2;

  double n_m0 = m0 + delta_t / 86400.0;
  double n_m1 = m1 + delta_t / 86400.0;
  double n_m2 = m2 + delta_t / 86400.0;

  double RA_interp_m0 = interpolate(RA_1, RA_2, RA_3, n_m0);
  double RA_interp_m1 = interpolate(RA_1, RA_2, RA_3, n_m1);
  double RA_interp_m2 = interpolate(RA_1, RA_2, RA_3, n_m2);

  double Decl_interp_m1 = interpolate(Decl_1, Decl_2, Decl_3, n_m1);
  double Decl_interp_m2 = interpolate(Decl_1, Decl_2, Decl_3, n_m2);

  double H_m0 = theta_0_m0 - LONGITUDE - RA_interp_m0;
  double H_m1 = theta_0_m1 - LONGITUDE - RA_interp_m1;
  double H_m2 = theta_0_m2 - LONGITUDE - RA_interp_m2;

  double h_m1 = std::asin(std::sin(LATITUDE * RAD) * std::sin(Decl_interp_m1 * RAD) + std::cos(LATITUDE * RAD) * std::cos(Decl_interp_m1 * RAD) * std::cos(H_m1 * RAD)) * DEG;
  double h_m2 = std::asin(std::sin(LATITUDE * RAD) * std::sin(Decl_interp_m2 * RAD) + std::cos(LATITUDE * RAD) * std::cos(Decl_interp_m2 * RAD) * std::cos(H_m2 * RAD)) * DEG;

  // transit correction
  double del_m0 = -1 * (H_m0 / 360);
  m0 += del_m0;

  // rise set correction
  double c1_m1 = std::cos(Decl_interp_m1 * RAD);
  double c2_m1 = std::cos(LATITUDE * RAD);
  double s1_m1 = std::sin(H_m1 * RAD);

  double c1_m2 = std::cos(Decl_interp_m2 * RAD);
  double c2_m2 = std::cos(LATITUDE * RAD);
  double s1_m2 = std::sin(H_m2 * RAD);

  double del_m1 = ((h_m1 - h_0) / (360 * c1_m1 * c2_m1 * s1_m1));
  double del_m2 = ((h_m2 - h_0) / (360 * c1_m2 * c2_m2 * s1_m2));

  m1 += del_m1;
  m2 += del_m2;

  m0 = constrain1(m0);
  m1 = constrain1(m1);
  m2 = constrain1(m2);

  std::string m0_s = timeDecToString(m0);
  std::string m1_s = timeDecToString(m1);
  std::string m2_s = timeDecToString(m2);
  if (DEBUG_RISE_SET)
  {
    printf("Transit: %lf, Rise: %lf, Set: %lf\n\n", m0, m1, m2);
    printf("h_m1: %lf, h_m2: %lf\n", h_m1, h_m2);
    printf("H_m0: %lf, H_m1: %lf, H_m2: %lf\n", H_m0, H_m1, H_m2);
    printf("del_m0: %lf, del_m1: %lf, del_m2: %lf\n", del_m0, del_m1, del_m2);
    printf("RA_1: %lf, RA_2: %lf, RA_3: %lf\n", RA_1, RA_2, RA_3);
    printf("Decl_1: %lf, Decl_2: %lf, Decl_3: %lf\n", Decl_1, Decl_2, Decl_3);
    printf("RA_interp_m0: %lf, RA_interp_m1: %lf, RA_interp_m2: %lf\n", RA_interp_m0, RA_interp_m1, RA_interp_m2);
    printf("Decl_interp_m1: %lf, Decl_interp_m2: %lf\n", Decl_interp_m1, Decl_interp_m2);
  }

  result.transit = m0_s;
  result.rise = m1_s;
  result.set = m2_s;

  return result;
}

Moon getMoonData()
{
  JDTime jd_time = julianDate();
  JDTime jd_time_p1 = julianDate(1);
  JDTime jd_time_m1 = julianDate(-1);
  Moon data;
  MoonPos pos;
  MoonPos pos_0;
  MoonPos pos_p1;
  MoonPos pos_m1;
  pos = getMoonPos(jd_time, 0);
  pos_0 = getMoonPos(jd_time, 1);
  pos_p1 = getMoonPos(jd_time_p1, 1);
  pos_m1 = getMoonPos(jd_time_m1, 1);

  // set past and present variables for interpolations

  pos_0.RA_p1 = pos_p1.RA;
  pos_0.Decl_p1 = pos_p1.Decl;
  pos_0.RA_m1 = pos_m1.RA;
  pos_0.Decl_m1 = pos_m1.Decl;

  MoonRiseSet rise_set;
  rise_set = getMoonRiseSet(pos_0, jd_time);

  const double i = pos.i;
  double d = moonPercent(i);
  std::string phase = getPhase(i, d);

  data.perc = d * 100;
  data.phase = phase;
  data.altitude = pos.altitude;
  data.azimuth = pos.azimuth;
  data.transit = rise_set.transit;
  data.rise = rise_set.rise;
  data.set = rise_set.set;
  return data;
}

int main()
{
  Moon m = getMoonData();
  printf("Transit: %s, Rise: %s, Set: %s\n", m.transit.c_str(), m.rise.c_str(), m.set.c_str());
  printf("Altitude: %lf\nAzimuth: %lf\n", m.altitude, m.azimuth);
  printf("%lf%%\n%s\n", m.perc, m.phase.c_str());
  return 0;
}
