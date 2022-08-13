#ifndef TABLE-CONSTANTS_H_INCLUDED
#define TABLE-CONSTANTS_H_INCLUDED

//lengths are in inches.

const int framerate=100;
const double timestep=1./framerate;

const double epsilon=pow(10.,-12.);

const int updates=floor(1./(timestep*framerate));

const double in2m=1./39.37;
const double pi=3.1415926535897932384626433832795028841971;
const double gravity=9.81/in2m;

const double table_length=140.25;
const double table_width=70.125;
const double cush_thickness=1.875;
const double rail_thickness=2.445;
const double bounciness=0.7;
const double mus=0.0075; //0.005-0.015
const double muk=0.15; //0.15-0.4

const double panel_ratio=0.2;
const double raw_width=(2.*rail_thickness+2.*cush_thickness+table_length);
const double raw_height=(2.*rail_thickness+2.*cush_thickness+table_width)*(1.+panel_ratio);
const double dfactor=sf::VideoMode::getDesktopMode().height/raw_height;
const double window_width=(2.*rail_thickness+2.*cush_thickness+table_length)*dfactor;
const double window_height=(2.*rail_thickness+2.*cush_thickness+table_width)*dfactor;
const double panel_height=window_height*panel_ratio;

const double ball_radius=0.02625/in2m;
const double ball_mass=0.141;
const double eball=1.; //coefficient of restitution.

const double mpockets[2][2]={{rail_thickness+cush_thickness+table_width,rail_thickness+2.*cush_thickness+table_width+0.156},
                             {rail_thickness+cush_thickness+table_width,rail_thickness-0.156}};
const double cpockets[4][2]={{rail_thickness,rail_thickness},
                             {rail_thickness,rail_thickness+2.*cush_thickness+table_width},
                             {rail_thickness+2.*cush_thickness+table_length,rail_thickness},
                             {rail_thickness+2.*cush_thickness+table_length,rail_thickness+2.*cush_thickness+table_width}};

const double spot_r=0.25;

const double mpocket_r=2.094;
const double cpocket_r=3.5;

const double pround=3.1;
const double k1=8.595-9.*sin(pi/4.);
const double k2=k1+4.5*sin(pi/4.)-2.445;

const double cueball_break_x=rail_thickness+cush_thickness+table_length-29.;
const double cueball_break_y=rail_thickness+cush_thickness+table_width/2.+5.3;

const double yellow_x=rail_thickness+cush_thickness+table_length-29.;
const double yellow_y=rail_thickness+cush_thickness+table_width/2.+11.687;
const double green_x=yellow_x;
const double green_y=rail_thickness+cush_thickness+table_width/2.-11.687;
const double brown_x=yellow_x;
const double brown_y=rail_thickness+cush_thickness+table_width/2.;
const double blue_x=rail_thickness+cush_thickness+table_width;
const double blue_y=brown_y;
const double pink_x=rail_thickness+cush_thickness+table_width/2.;
const double pink_y=brown_y;
const double black_x=rail_thickness+cush_thickness+12.75;
const double black_y=brown_y;

const double colourpos[6][2]={{yellow_x,yellow_y},{green_x,green_y},{brown_x,brown_y},{blue_x,blue_y},{pink_x,pink_y},{black_x,black_y}};

const sf::Color yellow_col=sf::Color(255,255,0);
const sf::Color green_col=sf::Color(0,150,0);
const sf::Color brown_col=sf::Color(131,87,43);
const sf::Color blue_col=sf::Color(0,0,255);
const sf::Color pink_col=sf::Color(255,105,180);
const sf::Color black_col=sf::Color(0,0,0);

//colours.
const int railcolour[3]={55,18,0};
const int baizecolour[3]={0,110,0};
const int cushioncolour[3]={0,80,0};
const int leathercolour[3]={255,229,153};

//cushion height-radius.
const double cushion_diff=0.03/in2m-ball_radius;
const double cushion_alpha=asin(cushion_diff/ball_radius);
const double cush_z1=0.03/in2m;
const double cush_z2=0.039/in2m;

//pocket things.
const double cpocket_angle=atan2(3.33402,1.06503);
const double mpocket_angle=atan2(1.93152,0.652752+0.156);

std::array<double,2> xlim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_length-ball_radius};
std::array<double,2> ylim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_width-ball_radius};

std::array<std::array<double,5>,3> cp={{{0.5*k1,0.5*k1,sqrt(2.)*0.5*k1-ball_radius,0.75*pi,-0.75*pi},
                                        {6.15,-2.445,4.5+ball_radius,-0.25*pi,-atan(6./11.)},
                                        {5.43,-1.125,3.+ball_radius,-atan(6./11.),0.}}};
std::array<std::array<double,5>,3> mp={{{0.,1.875,pround-ball_radius,pi-asin(1.625/pround),pi},
                                        {3.719,-0.438,2.094+ball_radius,-0.5*pi,-atan(0.625/0.812)},
                                        {4.344,-1.25,3.125+ball_radius,-atan(0.625/0.812),0.}}};

std::array<double,2> mpline={blue_x-1.625,blue_x+1.625};
std::array<std::array<double,4>,8> cpline={{{1.,-k1,rail_thickness+k1,rail_thickness+k2},
                                            {1.,k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+table_width-k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+table_width+k1,rail_thickness+k1,rail_thickness+k2},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+k1+table_length,rail_thickness+2.*cush_thickness+table_length-k2+k1,rail_thickness+2.*cush_thickness+table_length},
                                            {-1.,2.*rail_thickness+2.*cush_thickness-k1+table_length,rail_thickness+2.*cush_thickness+table_length-k2,rail_thickness+2.*cush_thickness+table_length-k1},
                                            {1.,-k1-table_width,rail_thickness+2.*cush_thickness+table_length-k2+k1,rail_thickness+2.*cush_thickness+table_length},
                                            {1.,k1-table_width,rail_thickness+2.*cush_thickness+table_length-k2,rail_thickness+2.*cush_thickness+table_length-k1}}};


#endif // TABLE-CONSTANTS_H_INCLUDED
