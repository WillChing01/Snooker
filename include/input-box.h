#ifndef INPUT-BOX_H_INCLUDED
#define INPUT-BOX_H_INCLUDED

const sf::Color inputOutlineColour1=sf::Color(200,200,200,255);
const sf::Color inputOutlineColour2=sf::Color(255,255,255,255);
const sf::Color inputColour1=sf::Color(100,100,100,150);

class InputBox
{
    public:
        double _textfactor=0.75;
        double _absoutlinethickness=2.;
        double _abscursorthickness=1.;

        double _t=0.;
        double _ton=0.8;
        double _toff=0.8;

        bool _isactive=true;
        bool _isTyping=false;

        int _cursorpos=0;

        std::string _input;

        sf::RectangleShape _shape;
        sf::RectangleShape _cursor;
        sf::Font _font;
        sf::Text _text;

        sf::Text _backtext;

        sf::Color _outlinecolour1=inputOutlineColour1;
        sf::Color _outlinecolour2=inputOutlineColour2;

        sf::Color _colour1=inputColour1;

        InputBox() {};
};

#endif // INPUT-BOX_H_INCLUDED
