//
// Created by Filippo Casari on 26.04.23.
//

#ifndef AS3PM_PARTICLE_H
#define AS3PM_PARTICLE_H


#include <cmath>
#include <iostream>
class Particle
{
public:
    double x;
    double y;
    double thermostat_temp;
    double vx;
    double vy;
    double m;
    double L;
    double a_x = 0;
    double a_y = 0;
    double F_x=0 ;
    double F_y=0 ;
    int i;
    int j;
    double k_en=0;
    Particle(double x, double y, double vx, double vy, double m, double L, int i, int j, double thermostat_temp)
    {
        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
        this->m = m;
        this->L = L;
        this->i = i;
        this->j = j;
        this->thermostat_temp = thermostat_temp;
        this->compute_k_energy();
    }
    void set_forces_zero(){
        this->F_x =0;
        this->F_y =0;
    }
    void update(double dt, double rc)
    {
        this->y += (this->vy * dt) + this->a_y * dt * dt * 0.5;
        this->x += (this->vx * dt) + this->a_x * dt * dt * 0.5;
        if(this->x < 0){
            this->x += this->L;
        }else if(this->x > this->L){
            this->x -= this->L;
        }
        if(this->y < 0){
            this->y += this->L;
        }else if(this->y > this->L){
            this->y -= this->L;
        }

        this->i = int(this->x / rc);
        this->j = int(this->y / rc);
        if(this->i<0){
            this->i += (this->L/rc);
        }
        if(this->j<0){
            this->j += (this->L/rc);
        }
    }
    void update_vel_acc(double dt, double Fx, double Fy)
    {

        this->F_x = Fx;
        this->F_y = Fy;
        double tmp_vel = sqrt(this->vx*this->vx + this->vy*this->vy);
        double new_acc_x = this->F_x / this->m;
        double new_acc_y = this->F_y / this->m;
        this->vx += (this->a_x + new_acc_x) * dt * 0.5;
        this->vy += (this->a_y + new_acc_y) * dt * 0.5;
        this->a_x = new_acc_x;
        this->a_y = new_acc_y;
        double new_vel = sqrt(this->vx*this->vx + this->vy*this->vy);
        //std::cout<<"difference velocity : "<< new_vel-tmp_vel<<std::endl;
        this->compute_k_energy();
    }
    void compute_k_energy()
    {
        this->k_en = 0.5 * this->m * (pow(this->vx, 2) + pow(this->vy, 2));
    }
    bool operator==(const Particle& other) const {
        return this->x == other.x && this->y == other.y && this->vx == other.vx && this->vy == other.vy && this->m == other.m && this->L == other.L && this->i == other.i && this->j == other.j;
    }
    bool operator!=(const Particle& other) const {
        return !(*this == other);
    }
    void apply_thermostat_scaling(double temp){
        double gamma = sqrt(1 + 0.0025*((this->thermostat_temp/temp)-1));
        this->vx *= gamma;
        this->vy *= gamma;
    }

};


#endif //AS3PM_PARTICLE_H
