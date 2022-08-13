#ifndef CUE_H_INCLUDED
#define CUE_H_INCLUDED

//cue texture file.
std::string cuetexturefile="cue0.png";

class Cue
{
    public:
        double _angle=-0.5*pi;
        double _speed=30.;
        double _mass=0.525;
        double _alpha=0.;
        double _offset=0.;
        double _theta=0.;
        double _eta=0.87;

        double _ballv;
        double _ballvz;
        double _ballparspin;
        double _ballperspin;
        double _ballrspin;

        boost::random_device rd;
        boost::mt19937 _gen{boost::mt19937(rd())};
        boost::normal_distribution<> _ndangle{boost::normal_distribution<>(0.0,pi/3600.)};
        boost::normal_distribution<> _ndalpha{boost::normal_distribution<>(0.0,0.003)};
        boost::normal_distribution<> _ndspeed{boost::normal_distribution<>(0.0,0.0085)};
        boost::normal_distribution<> _ndoffset{boost::normal_distribution<>(0.0,0.04)};

        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varangle{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndangle)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varalpha{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndalpha)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varspeed{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndspeed)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varoffset{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndoffset)};

        sf::Texture _texture;
        sf::Sprite _sprite;

        std::string tfile=cuetexturefile;

        Cue();
        void shot();
        void perturb();
};

Cue::Cue()
{
    if (!_texture.loadFromFile(tfile))
    {
        std::cout << "Error loading cue texture!" << std::endl;
    }
    _texture.setSmooth(true);
    _sprite.setTexture(_texture);
//    _sprite.scale(58.*dfactor/1369.,58.*dfactor/1369.);
//    _sprite.setOrigin(0.,26.*58.*dfactor/1369.);
    _sprite.scale(57.*dfactor/5213.,57.*dfactor/5213.);
    sf::FloatRect bounds=_sprite.getLocalBounds();
    _sprite.setOrigin(bounds.left,bounds.top+bounds.height*0.5);
}

void Cue::shot()
{
    double dy=_offset*cos(_theta);
    double dx=_offset*sin(_theta);
    double a=(1.+(ball_mass/_mass)+(2.5/pow(ball_radius,2.))*(pow(_offset,2.)));
    double v=(_speed/a)*(1+sqrt(_eta-(1-_eta)*(_mass/ball_mass)*(1+(2.5/pow(ball_radius,2.))*(pow(_offset,2.)))));
    double w=2.5*v*_offset/pow(ball_radius,2.);

    _ballparspin=w*cos(_theta);
    _ballrspin=-w*sin(_theta)*cos(_alpha);
    _ballperspin=w*sin(_theta)*sin(_alpha);

    _ballv=v*cos(_alpha);
    _ballvz=-v*sin(_alpha);
}

void Cue::perturb()
{
    try
    {
    _angle+=_varangle();
    _alpha+=_varalpha();
    _speed+=_varspeed();

    if (_alpha<0.) {_alpha=0.;}
    if (_alpha>0.5*pi) {_alpha=0.5*pi;}

    if (_speed<0) {_speed=0.;}
    if (_speed>100.5*1.2) {_speed=100.5*1.2;}

    double dx=_offset*sin(_theta);
    double dy=_offset*cos(_theta);

    dx+=_varoffset();
    dy+=_varoffset();

    _offset=sqrt(pow(dx,2.)+pow(dy,2.));
    _theta=atan2(dx,dy);
    }
    catch (...) {}

    if (_offset>ball_radius) {_offset=ball_radius;}
}

#endif // CUE_H_INCLUDED
