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
            creator.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            creator.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.2));
            creator.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&creator);

            bigtitle.setFont(_boldfont);
            bigtitle.setCharacterSize(int(_sfac*raw_height*0.15));
            bigtitle.setString("Snooker");
            textrect=bigtitle.getLocalBounds();
            bigtitle.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            bigtitle.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.3));
            bigtitle.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&bigtitle);

            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());

            double buttonwidth=_sfac*raw_width*0.15;
            int textcolour1[4]={255,255,255,255};
            int textcolour2[4]={255,255,255,255};
            int textcolour3[4]={255,255,255,255};
            int outlinecolour1[4]={200,200,200,255};
            int outlinecolour2[4]={255,255,255,255};
            int outlinecolour3[4]={255,255,255,255};
            int colour1[4]={100,100,100,150};
            int colour2[4]={100,100,100,150};
            int colour3[4]={169,169,169,200};

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
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));

                _buttons[i]._colour1=sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]);
                _buttons[i]._colour2=sf::Color(colour2[0],colour2[1],colour2[2],colour2[3]);
                _buttons[i]._colour3=sf::Color(colour3[0],colour3[1],colour3[2],colour3[3]);

                _buttons[i]._textcolour1=sf::Color(textcolour1[0],textcolour1[1],textcolour1[2],textcolour1[3]);
                _buttons[i]._textcolour2=sf::Color(textcolour2[0],textcolour2[1],textcolour2[2],textcolour2[3]);
                _buttons[i]._textcolour3=sf::Color(textcolour3[0],textcolour3[1],textcolour3[2],textcolour3[3]);

                _buttons[i]._outlinecolour1=sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]);
                _buttons[i]._outlinecolour2=sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]);
                _buttons[i]._outlinecolour3=sf::Color(outlinecolour3[0],outlinecolour3[1],outlinecolour3[2],outlinecolour3[3]);

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }
        }
        void update(double dt,sf::Vector2i mouse_pos) {};
};

#endif // TITLE-SCREEN_H_INCLUDED
