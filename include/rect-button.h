#ifndef RECT-BUTTON_H_INCLUDED
#define RECT-BUTTON_H_INCLUDED

class RectButton
{
    public:
        double _ratio=4.;
        //ratio is width/height.
        double _textfactor=0.75;
        double _absoutlinethickness=2.;

        bool _controlchange=false;
        bool _isactive=true;

        std::string _target;

        sf::RectangleShape _shape;
        sf::Font _font;
        sf::Text _text;

        sf::Color _textcolour1;
        sf::Color _textcolour2;
        sf::Color _textcolour3;

        sf::Color _colour1;
        sf::Color _colour2;
        sf::Color _colour3;

        sf::Color _outlinecolour1;
        sf::Color _outlinecolour2;
        sf::Color _outlinecolour3;

        RectButton() {};
};

#endif // RECT-BUTTON_H_INCLUDED
