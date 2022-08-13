#ifndef INPUT-BOX_H_INCLUDED
#define INPUT-BOX_H_INCLUDED

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

        int _cursorpos=0;

        std::string _input;

        sf::RectangleShape _shape;
        sf::RectangleShape _cursor;
        sf::Font _font;
        sf::Text _text;

        sf::Text _backtext;

        sf::Color _outlinecolour1;
        sf::Color _outlinecolour2;

        InputBox() {};
};

#endif // INPUT-BOX_H_INCLUDED
