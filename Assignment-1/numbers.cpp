// Copyright year Negmetulla_Yerlan 180107219@stu.sdu.edu.kz

#include <vector>
#include <string>
using std::vector;
using std::string;

bool is_happy(int x) {
  vector<int> x_set;
  bool flag = true;
  while (x > 1 && flag) {
    x_set.push_back(x);
    int total = 0;
    while (x > 0) {
      total += (x%10) * (x%10); x /= 10;
    }
    x = total;
    for (int e: x_set)
      if (x == e) {flag = false; break;}
  }
  return x == 1;
}

double product_of_positives(vector<double> vect) {
  double total = 1;
  for (auto i = vect.begin(); i != vect.end(); ++i) 
    if (*i > 0)
      total *= *i;
  return total;
}

int product_of_positives(vector<int> vect) {
  int total = 1;
  for (auto i = vect.begin(); i != vect.end(); ++i) 
    if (*i > 0)
      total *= *i;
  return total;
}

vector<int> proper_divisors(int n) {
  vector<int> vect;
  for (int i=1; i < n/2+1; i++) 
    if (n % i == 0) vect.push_back(i);
  return vect;
}

string add(const string& num1, const string& num2) {
  string result;
  result.reserve( 1 + std::max(num1.length(), num2.length()) );
  int num1pos = num1.length();
  int num2pos = num2.length();
  int carry = 0;
  while( carry > 0 || num1pos > 0 || num2pos > 0 ) {
    if( num1pos > 0 ) carry += num1.at(--num1pos) - '0';
    if( num2pos > 0 ) carry += num2.at(--num2pos) - '0';
    result.push_back('0' + (carry%10));
    carry /= 10;
  }
  for (int i=0; i<result.length()/2; i++) {
    char temp = result.at(i);
    result.at(i) = result.at(result.length()-1-i);
    result.at(result.length()-1-i) = temp;
  }
  return result;
}
