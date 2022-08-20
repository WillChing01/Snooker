#ifndef TITLE-SCREEN_H_INCLUDED
#define TITLE-SCREEN_H_INCLUDED

class TitleScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text creator;
        sf::Text bigtitle;

    public:
        TitleScreen(double sfac=dfactor) : GameState(sfac)
        {
            _background=sf::Color(0,0,0);

            if (!_thinfont.loadFromFile(_thinFontFile))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile(_boldFontFile))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            creator.setFont(_thinfont);
            creator.setCharacterSize(int(_sfac*raw_height*0.04));
            creator.setString("William Ching's");
            textrect=creator.getLocalBounds();
            creator.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            creator.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.2)));
            creator.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&creator);

            bigtitle.setFont(_boldfont);
            bigtitle.setCharacterSize(int(_sfac*raw_height*0.15));
            bigtitle.setString("Snooker");
            textrect=bigtitle.getLocalBounds();
            bigtitle.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            bigtitle.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.3)));
            bigtitle.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&bigtitle);

            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());

            double buttonwidth=_sfac*raw_width*0.15;

            for (int i=0;i<4;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Singleplayer");
                    _buttons[i]._target="Singleplayer";
                }
                else if (i==1)
                {
                    _buttons[i]._text.setString("Multiplayer");
                    _buttons[i]._target="Multiplayer";
                }
                else if (i==2)
                {
                    _buttons[i]._text.setString("Options");
                    _buttons[i]._target="Options";
                }
                else if (i==3)
                {
                    _buttons[i]._text.setString("Quit");
                    _buttons[i]._target="Quit";
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
                _buttons[i]._text.setPosition(sf::Vector2f(int(0.5*_sfac*raw_width),int(0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio)));

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }
        }
        void update(double dt,sf::Vector2i mouse_pos) {};
};

#endif // TITLE-SCREEN_H_INCLUDED
