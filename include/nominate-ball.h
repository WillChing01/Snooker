#ifndef NOMINATE-BALL_H_INCLUDED
#define NOMINATE-BALL_H_INCLUDED

class NominateBall
{
    public:
        sf::CircleShape _shape;
        double _abslinethickness=2.;
        sf::Color _colour;

        NominateBall() {};
};

#endif // NOMINATE-BALL_H_INCLUDED
