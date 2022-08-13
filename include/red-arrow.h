#ifndef RED-ARROW_H_INCLUDED
#define RED-ARROW_H_INCLUDED

class RedArrow
{
    public:
        sf::CircleShape _head;
        sf::RectangleShape _tail;

        RedArrow();
};

RedArrow::RedArrow()
{
    _head.setRadius(0.75*ball_radius*dfactor);
    _head.setOrigin(dfactor*0.75*ball_radius,dfactor*0.75*ball_radius);
    _head.setPosition(sf::Vector2f(-1000.,-1000.));
    _head.setPointCount(3);
    _head.rotate(180.);
    _tail.setSize(sf::Vector2f(0.7*ball_radius*dfactor,2*ball_radius*dfactor));
    _tail.setOrigin(dfactor*0.35*ball_radius,dfactor*ball_radius);
    _tail.setPosition(sf::Vector2f(-1000.,-1000.));

    _head.setFillColor(sf::Color(255,0,0));
    _tail.setFillColor(sf::Color(255,0,0));
}

#endif // RED-ARROW_H_INCLUDED
