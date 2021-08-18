// Copyright year Negmetulla_Yerlan 180107219@stu.sdu.edu.kz
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <limits>
#include <math.h>
#include <sys/stat.h>
using namespace std;


/*
This file for check results. You need "my_outputs" folder
How to run the file?

> g++ testing.cpp -o testing
> ./testing

*/


class Vector3d {
public:
    double x, y, z;
    Vector3d(double x, double y, double z) {
        this -> x = x;
        this -> y = y;
        this -> z = z;
    }
    Vector3d() {
        this -> x = 0;
        this -> y = 0;
        this -> z = 0;
    }
    Vector3d operator+(Vector3d pos) 
        {return Vector3d(x + pos.x, y + pos.y, z + pos.z);}
    Vector3d operator-(Vector3d pos) 
        {return Vector3d(x - pos.x, y - pos.y, z - pos.z);}
    double operator*(Vector3d pos) 
        {return double(x*pos.x + y*pos.y + z*pos.z);}

    Vector3d operator*(double d) 
        {return Vector3d(x*d, y*d, z*d);}
};

class Sphere {
private:
      double   mass;
      double   radius;
      Vector3d position;
      Vector3d velocity;
      string   name;
      int      bounces;
      bool     cont;
public:
    // Constructor
    Sphere(double mass, double radius, Vector3d position, Vector3d velocity, string name) {
        this -> mass = mass;
        this -> radius = radius;
        this -> position = position;
        this -> velocity = velocity;
        this -> name = name;
        this -> bounces = 0;
        this -> cont = false;
    }
    Sphere(double radius, string name, bool cont) {
        this -> radius = radius;
        this -> name = name;
        this -> mass = numeric_limits<double>::infinity();
        this -> position = Vector3d();
        this -> velocity = Vector3d();
        this -> bounces = 0;
        this -> cont = cont;
    }
    // Get
    string   getName()     {return this -> name;}
    double   getMass()     {return this -> mass;}
    double   getRadius()   {return this -> radius;}
    Vector3d getPosition() {return this -> position;}
    Vector3d getVelocity() {return this -> velocity;}
    int      getBounces()  {return this -> bounces;}
    // Set and Update
    void setPosition(Vector3d position) {this -> position = position;}
    void setVelocity(Vector3d velocity) {this -> velocity = velocity;}
    void updatePosition(double time) {this->position = this->position + this->velocity * time;}
    void updateVelocity(Sphere s) {
        Vector3d _v = this->calculateElasticCollision(s);
        this->setVelocity(_v);
    }
    void increaseBounces() {bounces++;}
    // Other
    bool checkCollisionCourse(Sphere s) {
        return (this->position - s.position) * (this->velocity - s.velocity) < 0;
    }
    double calculateDiscriminant(double a, double b, double c) {
        double zero = 1e-12;
        double discriminant = b*b - 4*a*c;
        if (discriminant > 0) {
            double t1 = (-b - sqrt(discriminant))/(2*a);
            double t2 = (-b + sqrt(discriminant))/(2*a);
            if (t1 > zero && t2 > zero)
                return t1 < t2 && t1 > zero ? t1 : t2;
            else if (t1 >  zero && t2 <= zero) return t1;
            else if (t1 <= zero && t2 >  zero) return t2;
            else return numeric_limits<double>::infinity();
        }
        else if (discriminant == 0) {
            double t = (-b)/(2*a);
            return t > zero ? t : numeric_limits<double>::infinity();
        }
        else return numeric_limits<double>::infinity();
    }
    double calculateCollisionTime(Sphere s) {
        Vector3d v1_v2 = this->velocity - s.velocity;
        Vector3d p1_p2 = this->position - s.position;
        double g = !s.cont ? 
            (this->radius + s.radius) * (this->radius + s.radius) : 
            (s.radius - this->radius) * (s.radius - this->radius);
        double a = v1_v2*v1_v2;
        double b = v1_v2*p1_p2 * 2;
        double c = p1_p2*p1_p2 - g;
        return this->calculateDiscriminant(a, b, c);
    }
    Vector3d calculateElasticCollision(Sphere s) {
        Vector3d v1_v2 = this->velocity - s.velocity;
        Vector3d r1_r2 = this->position - s.position;
        double g = !s.cont ? s.mass/(this->mass + s.mass) : 1;
        Vector3d _v = this->velocity - r1_r2*((v1_v2*r1_r2)/(r1_r2*r1_r2)) * 2*g;
        return _v;
    }
};



void updateVelocities(Sphere &s1, Sphere &s2) {
    Vector3d _v1 = s1.calculateElasticCollision(s2);
    Vector3d _v2 = s2.calculateElasticCollision(s1);
    s1.setVelocity(_v1);
    s2.setVelocity(_v2);
}
vector<Sphere> read_file(string filename) {
    fstream file;
    file.open(filename.c_str());
    if (!file.is_open()) 
        throw std::invalid_argument("File not found.");

    vector<Sphere> spheres;
    bool flag = true;

    string line;
    while(getline(file, line)) {
        flag = false;
        string word;
        istringstream iss0(line);

        int r=0;
        while(iss0 >> word) r++;
        if (r!=9) {flag = true;break;}

        istringstream iss(line);

        iss >> word; double  mass;
        if (stod(word) <= 0) {flag = true;break;}
        else mass = stod(word);
        iss >> word; double  radius;
        if (stod(word) <= 0) {flag = true;break;}
        else radius = stod(word);
        iss >> word; double  pos_x = stod(word);
        iss >> word; double  pos_y = stod(word);
        iss >> word; double  pos_z = stod(word);
        iss >> word; double  vel_x = stod(word);
        iss >> word; double  vel_y = stod(word);
        iss >> word; double  vel_z = stod(word);
        iss >> word; string  name = word;
        spheres.push_back(Sphere(mass, radius, Vector3d(pos_x, pos_y, pos_z), Vector3d(vel_x, vel_y, vel_z), name));
    }
    if (flag) throw std::invalid_argument("The file does not match the input format.");
    file.close();
    return spheres;
}



int write_predict_movements_and_collisions(double universe_radius, int bounces, string input_filename, string output_filename) {
    double total_time = 0;
    double soonest_time = numeric_limits<double>::infinity();
    vector<int> reflected_spheres;
    vector<array<int, 2> > collided_spheres;

    Sphere universe(universe_radius, "the universe", true);
    vector<Sphere> spheres;

    try {
        spheres = read_file(input_filename);
        ofstream file(output_filename.c_str());
        file << "Here are the initial conditions." << endl;
        file << "universe radius " << universe_radius << (universe_radius == (int)universe_radius ? ".0" : "") << endl;
        file << "max collisions " << bounces << endl;
        for (Sphere s: spheres) {
            string name = s.getName();
            double mass = s.getMass();
            double radius = s.getRadius();
            Vector3d pos = s.getPosition();
            Vector3d vel = s.getVelocity();
            int bounces = s.getBounces();
            file << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
        }

        double enegry = 0;
        for (Sphere s: spheres)
            enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
        file << "energy: " << enegry/2 << endl;

        Vector3d momentum = Vector3d();
        for (Sphere s: spheres)
            momentum = momentum + s.getVelocity() * s.getMass();
        file << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z << ")" << endl << endl;
        
        file << "Here are the events." << endl;

        while (spheres.size() > 0) {
            soonest_time = numeric_limits<double>::infinity();

            for (int i=0; i<spheres.size(); i++) {
                double time = spheres[i].calculateCollisionTime(universe);
                if (soonest_time == time) {
                    reflected_spheres.push_back(i);
                } else if (soonest_time > time) {
                    reflected_spheres.clear();
                    collided_spheres.clear();
                    reflected_spheres.push_back(i);
                    soonest_time = time;
                }
            }

            for (int i=0; i<spheres.size(); i++) {
                int n = 0;
                for (int j=i+1; j<spheres.size(); j++) 
                    if (spheres[i].checkCollisionCourse(spheres[j])) {
                        double time = spheres[i].calculateCollisionTime(spheres[j]);
                        if (soonest_time == time) {
                            collided_spheres.push_back({i,j});
                            n++;
                        } else if (soonest_time > time) {
                            reflected_spheres.clear();
                            collided_spheres.clear();
                            collided_spheres.push_back({i,j});
                            soonest_time = time;
                            n=1;
                        }
                        if (n > 1) break;
                    }
            }

            if (soonest_time == numeric_limits<double>::infinity()) break;
            for (Sphere &s: spheres) s.updatePosition(soonest_time);
            total_time += soonest_time;

            for (int i=0; i<reflected_spheres.size(); i++) {
                file << endl << "time of event: " << total_time << endl;
                spheres[reflected_spheres[i]].updateVelocity(universe);
                spheres[reflected_spheres[i]].increaseBounces();
                file << "reflecting " << spheres[reflected_spheres[i]].getName() << endl;

                for (Sphere s: spheres) {
                    string name = s.getName();
                    double mass = s.getMass();
                    double radius = s.getRadius();
                    Vector3d pos = s.getPosition();
                    Vector3d vel = s.getVelocity();
                    int bounces = s.getBounces();
                    file << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
                }
                enegry = 0;
                for (Sphere s: spheres)
                    enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
                file << "energy: " << enegry/2 << endl;

                momentum = Vector3d();
                for (Sphere s: spheres)
                    momentum = momentum + s.getVelocity() * s.getMass();
                file << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z << ")" << endl;
            }

            for (int i=0; i<collided_spheres.size(); i++) {
                file << endl << "time of event: " << total_time << endl;
                updateVelocities(spheres[collided_spheres[i][0]], spheres[collided_spheres[i][1]]);
                spheres[collided_spheres[i][0]].increaseBounces();
                spheres[collided_spheres[i][1]].increaseBounces();
                file << "colliding " << spheres[collided_spheres[i][0]].getName() << " " << spheres[collided_spheres[i][1]].getName() << endl;

                for (Sphere s: spheres) {
                    string name = s.getName();
                    double mass = s.getMass();
                    double radius = s.getRadius();
                    Vector3d pos = s.getPosition();
                    Vector3d vel = s.getVelocity();
                    int bounces = s.getBounces();
                    file << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
                }
                enegry = 0;
                for (Sphere s: spheres)
                    enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
                file << "energy: " << enegry/2 << endl;

                momentum = Vector3d();
                for (Sphere s: spheres)
                    momentum = momentum + s.getVelocity() * s.getMass();
                file << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z<<")" << endl;
            }

            int idx = 0;
            vector<Sphere> tmp_spheres = spheres;
            for (auto& e : tmp_spheres) {
                if (e.getBounces() >= bounces) {
                    file << endl << "disappear " << e.getName() << endl;
                    spheres.erase(spheres.begin() + idx);
                    idx = idx - 1;
                }
                idx = idx + 1;
            }
        }
        return 0;
    }
    catch(const std::invalid_argument& exc) {
        cerr << exc.what() << endl;
        return 1;
    }
}

int print_predict_movements_and_collisions(double universe_radius, int bounces, string input_filename) {
    double total_time = 0;
    double soonest_time = numeric_limits<double>::infinity();
    vector<int> reflected_spheres;
    vector<array<int, 2> > collided_spheres;

    Sphere universe(universe_radius, "the universe", true);
    vector<Sphere> spheres;

    try {
        if (input_filename == "") {
            bool flag = true;
            string line;
            while(getline(cin, line)) {
                flag = false;
                string word;
                istringstream iss0(line);

                int r=0;
                while(iss0 >> word) r++;
                if (r!=9) {flag = true;break;}

                istringstream iss(line);

                iss >> word; double  mass;
                if (stod(word) <= 0) {flag = true;break;}
                else mass = stod(word);
                iss >> word; double  radius;
                if (stod(word) <= 0) {flag = true;break;}
                else radius = stod(word);
                iss >> word; double  pos_x = stod(word);
                iss >> word; double  pos_y = stod(word);
                iss >> word; double  pos_z = stod(word);
                iss >> word; double  vel_x = stod(word);
                iss >> word; double  vel_y = stod(word);
                iss >> word; double  vel_z = stod(word);
                iss >> word; string  name = word;
                spheres.push_back(Sphere(mass, radius, Vector3d(pos_x, pos_y, pos_z), Vector3d(vel_x, vel_y, vel_z), name));
            }
            if (flag) throw std::invalid_argument("The input does not match the input format.");
        } else {
            spheres = read_file(input_filename);
        }
        cout << "Here are the initial conditions." << endl;
        cout << "universe radius " << universe_radius << (universe_radius == (int)universe_radius ? ".0" : "") << endl;
        cout << "max collisions " << bounces << endl;
        for (Sphere s: spheres) {
            string name = s.getName();
            double mass = s.getMass();
            double radius = s.getRadius();
            Vector3d pos = s.getPosition();
            Vector3d vel = s.getVelocity();
            int bounces = s.getBounces();
            cout << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
        }

        double enegry = 0;
        for (Sphere s: spheres)
            enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
        cout << "energy: " << enegry/2 << endl;

        Vector3d momentum = Vector3d();
        for (Sphere s: spheres)
            momentum = momentum + s.getVelocity() * s.getMass();
        cout << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z << ")" << endl << endl;
        
        cout << "Here are the events." << endl;

        while (spheres.size() > 0) {
            soonest_time = numeric_limits<double>::infinity();

            for (int i=0; i<spheres.size(); i++) {
                double time = spheres[i].calculateCollisionTime(universe);
                if (soonest_time == time) {
                    reflected_spheres.push_back(i);
                } else if (soonest_time > time) {
                    reflected_spheres.clear();
                    collided_spheres.clear();
                    reflected_spheres.push_back(i);
                    soonest_time = time;
                }
            }

            for (int i=0; i<spheres.size(); i++) {
                int n = 0;
                for (int j=i+1; j<spheres.size(); j++) 
                    if (spheres[i].checkCollisionCourse(spheres[j])) {
                        double time = spheres[i].calculateCollisionTime(spheres[j]);
                        if (soonest_time == time) {
                            collided_spheres.push_back({i,j});
                            n++;
                        } else if (soonest_time > time) {
                            reflected_spheres.clear();
                            collided_spheres.clear();
                            collided_spheres.push_back({i,j});
                            soonest_time = time;
                            n=1;
                        }
                        if (n > 1) break;
                    }
            }

            if (soonest_time == numeric_limits<double>::infinity()) break;
            for (Sphere &s: spheres) s.updatePosition(soonest_time);
            total_time += soonest_time;

            for (int i=0; i<reflected_spheres.size(); i++) {
                cout << endl << "time of event: " << total_time << endl;
                spheres[reflected_spheres[i]].updateVelocity(universe);
                spheres[reflected_spheres[i]].increaseBounces();
                cout << "reflecting " << spheres[reflected_spheres[i]].getName() << endl;

                for (Sphere s: spheres) {
                    string name = s.getName();
                    double mass = s.getMass();
                    double radius = s.getRadius();
                    Vector3d pos = s.getPosition();
                    Vector3d vel = s.getVelocity();
                    int bounces = s.getBounces();
                    cout << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
                }
                enegry = 0;
                for (Sphere s: spheres)
                    enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
                cout << "energy: " << enegry/2 << endl;

                momentum = Vector3d();
                for (Sphere s: spheres)
                    momentum = momentum + s.getVelocity() * s.getMass();
                cout << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z << ")" << endl;
            }

            for (int i=0; i<collided_spheres.size(); i++) {
                cout << endl << "time of event: " << total_time << endl;
                updateVelocities(spheres[collided_spheres[i][0]], spheres[collided_spheres[i][1]]);
                spheres[collided_spheres[i][0]].increaseBounces();
                spheres[collided_spheres[i][1]].increaseBounces();
                cout << "colliding " << spheres[collided_spheres[i][0]].getName() << " " << spheres[collided_spheres[i][1]].getName() << endl;

                for (Sphere s: spheres) {
                    string name = s.getName();
                    double mass = s.getMass();
                    double radius = s.getRadius();
                    Vector3d pos = s.getPosition();
                    Vector3d vel = s.getVelocity();
                    int bounces = s.getBounces();
                    cout << name<<" m="<<mass<<" R="<<radius<<" p=("<<pos.x<<","<<pos.y<<","<<pos.z<<") v=("<<vel.x<<","<<vel.y<<","<<vel.z<<") bounces="<<bounces<<endl;
                }
                enegry = 0;
                for (Sphere s: spheres)
                    enegry += (s.getVelocity()*s.getVelocity()) * s.getMass();
                cout << "energy: " << enegry/2 << endl;

                momentum = Vector3d();
                for (Sphere s: spheres)
                    momentum = momentum + s.getVelocity() * s.getMass();
                cout << "momentum: ("<< momentum.x<<","<<momentum.y<<","<<momentum.z<<")" << endl;
            }

            int idx = 0;
            vector<Sphere> tmp_spheres = spheres;
            for (auto& e : tmp_spheres) {
                if (e.getBounces() >= bounces) {
                    cout << endl << "disappear " << e.getName() << endl;
                    spheres.erase(spheres.begin() + idx);
                    idx = idx - 1;
                }
                idx = idx + 1;
            }
        }
        return 0;
    }
    catch(const std::invalid_argument& exc) {
        cerr << exc.what() << endl;
        return 1;
    }
}


bool compareFiles(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);
    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }
    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }
    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
}


int main(int argc, char const *argv[]) {
    if (mkdir("my_outputs", 0777) != -1)
        cout << "Directory created" << endl;

    struct sample{
        double universe_radius;
        int bounces;
        string input_filename;
        string output_filename;
        string check_filename;
    };
    vector<sample> vect;
    vect.push_back(sample{   100, 2, "examples_spheres/diagonal.txt",     "my_outputs/diagonal_100_2_JBC.txt",     "examples_spheres/diagonal_100_2_JBC.txt"     });
    vect.push_back(sample{    40, 2, "examples_spheres/glancing.txt",     "my_outputs/glancing_40_2_JBC.txt",      "examples_spheres/glancing_40_2_JBC.txt"      });
    vect.push_back(sample{  2000, 2, "examples_spheres/lineup.txt",       "my_outputs/lineup_2000_2_JBC.txt",      "examples_spheres/lineup_2000_2_JBC.txt"      });
    vect.push_back(sample{   200, 5, "examples_spheres/randompool.txt",   "my_outputs/randompool_200_5_JBC.txt",   "examples_spheres/randompool_200_5_JBC.txt"   });
    vect.push_back(sample{  2000, 2, "examples_spheres/slow.txt",         "my_outputs/slow_2000_2_JBC.txt",        "examples_spheres/slow_2000_2_JBC.txt"        });
    vect.push_back(sample{  2000, 5, "examples_spheres/smashup.txt",      "my_outputs/smashup_2000_5_JBC.txt",     "examples_spheres/smashup_2000_5_JBC.txt"     });
    vect.push_back(sample{ 10.41, 2, "examples_spheres/tennisball.txt",   "my_outputs/tennisball_10.41_2_JBC.txt", "examples_spheres/tennisball_10.41_2_JBC.txt" });
    vect.push_back(sample{   100, 3, "examples_spheres/tennisball.txt",   "my_outputs/tennisball_100_3_JBC.txt",   "examples_spheres/tennisball_100_3_JBC.txt"   });
    vect.push_back(sample{   120, 2, "examples_spheres/three.txt",        "my_outputs/three_120_2_JBC.txt",        "examples_spheres/three_120_2_JBC.txt"        });
    vect.push_back(sample{    50, 4, "examples_spheres/threeonxaxis.txt", "my_outputs/threeonxaxis_50_4_JBC.txt",  "examples_spheres/threeonxaxis_50_4_JBC.txt"  });
    vect.push_back(sample{   200, 2, "examples_spheres/xyzmovers.txt",    "my_outputs/xyzmovers_200_2_JBC.txt",    "examples_spheres/xyzmovers_200_2_JBC.txt"    });
    vect.push_back(sample{   100, 3, "examples_spheres/zmover.txt",       "my_outputs/zmover_100_3_JBC.txt",       "examples_spheres/zmover_100_3_JBC.txt"       });

    for (int i=0; i<vect.size(); i++) {
        write_predict_movements_and_collisions(vect[i].universe_radius, vect[i].bounces, vect[i].input_filename, vect[i].output_filename);
        cout << (compareFiles(vect[i].output_filename, vect[i].check_filename) ? "Correct -> " + vect[i].output_filename : "Incorrect -> " + vect[i].input_filename) << endl;
    }
    return 0;
}