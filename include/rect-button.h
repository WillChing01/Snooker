#ifndef RECT-BUTTON_H_INCLUDED
#define RECT-BUTTON_H_INCLUDED

const sf::Color buttonTextColour1=sf::Color(255,255,255,255);
const sf::Color buttonTextColour2=sf::Color(255,255,255,255);
const sf::Color buttonTextColour3=sf::Color(255,255,255,255);
const sf::Color buttonOutlineColour1=sf::Color(200,200,200,255);
const sf::Color buttonOutlineColour2=sf::Color(255,255,255,255);
const sf::Color buttonOutlineColour3=sf::Color(255,255,255,255);
const sf::Color buttonColour1=sf::Color(100,100,100,150);
const sf::Color buttonColour2=sf::Color(100,100,100,150);
const sf::Color buttonColour3=sf::Color(169,169,169,200);

//colour 1: default, no action.
//colour 2: when hover with mouse.
//colour 3: when left-clicked with mouse.

//at the moment, buttonColour3 must be unique from buttonColour1 and buttonColour2 !!!.

class RectButton
{
    public:
        double _ratio=4.;
        //ratio is width/height.
        double _textfactor=0.75;
        double _absoutlinethickness=absOutlineThickness;

        bool _controlchange=false; //for changing controls in options.
        bool _isactive=true;
        bool _wasClicked=false;
        bool _isClicked=false;
        bool _shouldExecute=false;

        std::string _target;

        sf::RectangleShape _shape;
        sf::Font _font;
        sf::Text _text;

        sf::Color _textcolour1=buttonTextColour1;
        sf::Color _textcolour2=buttonTextColour2;
        sf::Color _textcolour3=buttonTextColour3;

        sf::Color _colour1=buttonColour1;
        sf::Color _colour2=buttonColour2;
        sf::Color _colour3=buttonColour3;

        sf::Color _outlinecolour1=buttonOutlineColour1;
        sf::Color _outlinecolour2=buttonOutlineColour2;
        sf::Color _outlinecolour3=buttonOutlineColour3;

        RectButton() {};
        RectButton(double width, sf::Vector2f pos, std::string text, std::string target, sf::Font font);
};

RectButton::RectButton(double buttonwidth, sf::Vector2f pos, std::string text, std::string target, sf::Font font)
{
    _shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_ratio));
    _shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_ratio);
    _shape.setOutlineThickness(_absoutlinethickness);
    _shape.setPosition(pos);

    _text.setFont(font);
    _text.setCharacterSize(int(buttonwidth*_textfactor/_ratio));

    _text.setString(text);
    _target=target;

    sf::FloatRect textrect=_text.getLocalBounds();
    _text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
    _text.setPosition(pos);

    _shape.setFillColor(_colour1);
    _shape.setOutlineColor(_outlinecolour1);
    _text.setFillColor(_textcolour1);
}

#endif // RECT-BUTTON_H_INCLUDED
