//////////////////////////////////////////////////////////////////////////////
//
//  --- Object.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////


#include "common.h"

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Sphere::intersect(vec4 p0_w, vec4 V_w)
{
  IntersectionValues result;
  //TODO: Ray-sphere setup
    vec4 p0_o=this->INVC*p0_w;//p0 is the R_start
    vec4 V_o= (this->INVCStar * V_w);//V is the R_direction
    double mag = length(V_o);
    V_o = normalize(V_o);
    result.t_o=raySphereIntersection(p0_o, V_o);
    result.t_w=result.t_o/mag;
    result.P_w=p0_w+result.t_w*V_w;
    result.P_o=p0_o+result.t_o*V_o;
    result.N_o=normalize(result.P_o-vec4(0, 0, 0, 1));
    result.N_w=TRANINVC*result.N_o;
    result.N_w.w = 0.0;
    result.N_w = normalize(result.N_w);
  return result;
}

/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin O and radius r    : Find t ------- */
double Sphere::raySphereIntersection(vec4 p0, vec4 V, vec4 O, double r)
{
  //TODO: Ray-sphere intersection;
    double a=1.0;
    double b= dot(2.0*V,p0-O);
    double c=length(p0-O)*length(p0-O)-r*r;
    double in_sqrt=(b*b-4.0*a*c);
    if(in_sqrt<0)
    {
        return std::numeric_limits<double>::infinity();
    }
    double t_plus=(-b+sqrt(in_sqrt))/(2*a);
    double t_minus=(-b-sqrt(in_sqrt))/(2*a);
    if (t_plus<EPSILON)
    {
        t_plus=std::numeric_limits<double>::infinity();
    }
    if(t_minus<EPSILON)
    {
        t_minus=std::numeric_limits<double>::infinity();
    }
    return fmin(t_plus,t_minus);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Square::intersect(vec4 p0_w, vec4 V_w)
{
  IntersectionValues result;
  //TODO: Ray-square setup
    vec4 p0_o=this->INVC*p0_w;//p is the R_start
    vec4 V_o= (this->INVCStar * V_w);//V is the R_direction
    double mag = length(V_o);
    V_o = normalize(V_o);
    result.t_o=raySquareIntersection(p0_o, V_o);
    result.t_w=result.t_o/mag;
    result.P_w=p0_w+result.t_w*V_w;
    result.P_o=p0_o+result.t_o*V_o;
    result.N_o=(vec4(0, 0, 1, 0));
    result.N_w=TRANINVC*result.N_o;
    result.N_w.w = 0.0;
    result.N_w = normalize(result.N_w);
  return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(vec4 p0, vec4 V)
{
    vec4 s= vec4(0,0,0,1);
    vec4 surf_normal= vec4(0,0,1,0);
    double check=dot(surf_normal,V);
    if(check==0)
    {
        return std::numeric_limits<double>:: infinity();
    }
    else{
        double t=dot(surf_normal,(s-p0))/check;
        vec4 checkPoint= p0+t*V;
        if(t<EPSILON)
        {
            return std::numeric_limits<double>:: infinity();
            
        }
        if(checkPoint.x > -1 && checkPoint.x<1 &&checkPoint.y> -1 && checkPoint.y<1)
        {
            return t;
            
        }
    }
    return std::numeric_limits<double>:: infinity();
}
