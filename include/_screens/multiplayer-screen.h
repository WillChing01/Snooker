#ifndef MULTIPLAYER-SCREEN_H_INCLUDED
#define MULTIPLAYER-SCREEN_H_INCLUDED

class MultiplayerScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;
        sf::Text ip;
        sf::Text port;

    public:
        MultiplayerScreen(double sfac=dfactor,std::string ipad="N/A",std::string portnum="N/A") : GameState(sfac)
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

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Multiplayer");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            title.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.15)));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

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

            ip.setFont(_thinfont);
            ip.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
            ip.setString("IP : "+ipad);
            textrect=ip.getLocalBounds();
            ip.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            ip.setPosition(sf::Vector2f(int(_sfac*raw_width*0.3),int(_sfac*raw_height*0.5)));
            ip.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&ip);

            port.setFont(_thinfont);
            port.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
            port.setString("Port : "+portnum);
            textrect=port.getLocalBounds();
            port.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            port.setPosition(sf::Vector2f(int(_sfac*raw_width*0.3),int(_sfac*raw_height*0.5+1.2*buttonwidth/_buttons[0]._ratio)));
            port.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&port);

            for (int i=0;i<3;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Host game");
                    _buttons[i]._target="MultiplayerHost";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.3*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.3*_sfac*raw_width),int(0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio)));
                }
                else if (i==1)
                {
                    _buttons[i]._text.setString("Join game");
                    _buttons[i]._target="MultiplayerJoin";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.7*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.7*_sfac*raw_width),int(0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio)));
                }
                else if (i==2)
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.5*_sfac*raw_width),int(0.50*_sfac*raw_height+(4.*1.2)*buttonwidth/_buttons[i]._ratio)));
                }

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

            _inputboxes.push_back(InputBox());
            _inputboxes.push_back(InputBox());
            _inputboxes.push_back(InputBox());

            double boxwidth=_sfac*raw_width*0.3;

            for (int i=0;i<3;i++)
            {
                _inputboxes[i]._shape.setSize(sf::Vector2f(boxwidth,buttonwidth/_buttons[0]._ratio));
                _inputboxes[i]._shape.setOrigin(sf::Vector2f(0.5*boxwidth,0.5*buttonwidth/_buttons[0]._ratio));
                _inputboxes[i]._shape.setOutlineThickness(_inputboxes[i]._absoutlinethickness);
                _inputboxes[i]._shape.setPosition(sf::Vector2f(0.7*_sfac*raw_width,_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio));

                if (!_inputboxes[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _inputboxes[i]._text.setFont(_inputboxes[i]._font);
                _inputboxes[i]._text.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
                _inputboxes[i]._text.setFillColor(sf::Color(255,255,255));
                _inputboxes[i]._text.setString("");
                _inputboxes[i]._text.setPosition(sf::Vector2f(int(0.7*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness),int(_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio-0.5*int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio))));

                if (i==0)
                {
                    _inputboxes[i]._backtext.setString("Enter target IP...");
                }
                else if (i==1)
                {
                    _inputboxes[i]._backtext.setString("Enter target port...");
                }
                else if (i==2)
                {
                    _inputboxes[i]._backtext.setString("Enter name...");
                }
                _inputboxes[i]._backtext.setFont(_inputboxes[i]._font);
                _inputboxes[i]._backtext.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
                textrect=_inputboxes[i]._backtext.getLocalBounds();
                _inputboxes[i]._backtext.setOrigin(sf::Vector2f(int(textrect.left),int(textrect.top+textrect.height*0.5)));
                _inputboxes[i]._backtext.setFillColor(sf::Color(255,255,255,150));
                _inputboxes[i]._backtext.setPosition(sf::Vector2f(int(0.7*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness),int(_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio)));

                _inputboxes[i]._cursor.setSize(sf::Vector2f(_inputboxes[i]._abscursorthickness,buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
                _inputboxes[i]._cursor.setOrigin(sf::Vector2f(0.5*_inputboxes[i]._abscursorthickness,0.5*buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
                _inputboxes[i]._cursor.setPosition(sf::Vector2f(0.7*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio));

                _inputboxes[i]._outlinecolour1=sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]);
                _inputboxes[i]._outlinecolour2=sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]);

                _inputboxes[i]._shape.setFillColor(sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]));
                _inputboxes[i]._shape.setOutlineColor(_inputboxes[i]._outlinecolour1);
                _inputboxes[i]._cursor.setFillColor(sf::Color(255,255,255,0));

                if (i==2)
                {
                    _inputboxes[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio));
                    _inputboxes[i]._text.setPosition(sf::Vector2f(int(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness),int(_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio-0.5*int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio))));
                    _inputboxes[i]._backtext.setPosition(sf::Vector2f(int(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness),int(_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio)));
                    _inputboxes[i]._cursor.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio));
                }
            }
        }
        void update(double dt,sf::Vector2i mouse_pos)
        {
            for (int i=0;i<_inputboxes.size();i++)
            {
                if (_inputboxes[i]._shape.getOutlineColor()==_inputboxes[i]._outlinecolour2)
                {
                    _inputboxes[i]._t+=dt;
                    if (int(_inputboxes[i]._cursor.getFillColor().a)==255 && _inputboxes[i]._t>_inputboxes[i]._ton)
                    {
                        _inputboxes[i]._t=0.;
                        sf::Color c=_inputboxes[i]._cursor.getFillColor();
                        _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,0));
                    }
                    else if (int(_inputboxes[i]._cursor.getFillColor().a)==0 && _inputboxes[i]._t>_inputboxes[i]._toff)
                    {
                        _inputboxes[i]._t=0.;
                        sf::Color c=_inputboxes[i]._cursor.getFillColor();
                        _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));
                    }
                }
            }
        };
};

#endif // MULTIPLAYER-SCREEN_H_INCLUDED
