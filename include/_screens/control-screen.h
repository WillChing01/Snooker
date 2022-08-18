#ifndef CONTROL-SCREEN_H_INCLUDED
#define CONTROL-SCREEN_H_INCLUDED

class ControlScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;

        std::vector<sf::Text> t;

    public:
        ControlScreen(double sfac=dfactor) : GameState(sfac)
        {
            _shouldUpdate=false;

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
            title.setString("Controls");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            title.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.15)));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

            for (int i=0;i<default_controls.size();i++)
            {
                sf::Text example;
                t.push_back(example);
            }

            for (int i=0;i<default_controls.size()+2;i++)
            {
                _buttons.push_back(RectButton());
            }

            double buttonwidth=_sfac*raw_width*0.083;
            double height=0.35*_sfac*raw_height;
            int textcolour1[4]={255,255,255,255};
            int textcolour2[4]={255,255,255,255};
            int textcolour3[4]={255,255,255,255};
            int outlinecolour1[4]={200,200,200,255};
            int outlinecolour2[4]={255,255,255,255};
            int outlinecolour3[4]={255,255,255,255};
            int colour1[4]={100,100,100,150};
            int colour2[4]={100,100,100,150};
            int colour3[4]={169,169,169,200};

            const std::vector<std::string> controlsOrder=
            {
                "Pause game",
                "Toggle mute",
                "Aim left",
                "Aim right",
                "Precise aim left",
                "Precise aim right",
                "Increase power",
                "Decrease power",
                "Increase cue elevation",
                "Decrease cue elevation",
                "Offset cue tip up",
                "Offset cue tip down",
                "Offset cue tip left",
                "Offset cue tip right",
                "Strike cueball",
                "Move ball up",
                "Move ball down",
                "Move ball left",
                "Move ball right",
                "Place ball"
            };

            for (int i=0;i<default_controls.size();i++)
            {
                t[i].setFont(_thinfont);
                t[i].setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));

                _buttons[i]._controlchange=true;
                _buttons[i]._shape.setSize(sf::Vector2f(4.*buttonwidth/_buttons[i]._ratio,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(2.*buttonwidth/_buttons[i]._ratio,0.5*buttonwidth/_buttons[i]._ratio);

                if (i<default_controls.size()/2)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.475*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.775*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+((i-default_controls.size()/2)*1.2)*buttonwidth/_buttons[i]._ratio));
                }

                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));

                t[i].setString(controlsOrder[i]);
                _buttons[i]._target=controlsOrder[i];
                _buttons[i]._text.setString(KeyToString(user_controls[controlsOrder[i]]));

                textrect=t[i].getLocalBounds();
                t[i].setOrigin(sf::Vector2f(int(textrect.left),int(textrect.top)));
                t[i].setFillColor(getControlColour(t[i].getString()));

                textrect=_buttons[i]._text.getLocalBounds();

                _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));

                if (i<default_controls.size()/2)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.475*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio),int(height+(i*1.2)*buttonwidth/_buttons[i]._ratio)));
                    t[i].setPosition(sf::Vector2f(int(0.225*_sfac*raw_width),int(height+(i*1.2)*buttonwidth/_buttons[i]._ratio-0.5*buttonwidth/_buttons[i]._ratio+0.5*0.25*buttonwidth/_buttons[i]._ratio)));
                }
                else
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.775*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio),int(height+((i-default_controls.size()/2)*1.2)*buttonwidth/_buttons[i]._ratio)));
                    t[i].setPosition(sf::Vector2f(int(0.525*_sfac*raw_width),int(height+((i-default_controls.size()/2)*1.2)*buttonwidth/_buttons[i]._ratio-0.5*buttonwidth/_buttons[i]._ratio+0.5*0.25*buttonwidth/_buttons[i]._ratio)));
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

                _shapes.push_back(&t[i]);
            }

            buttonwidth=_sfac*raw_width*0.15;
            for (int i=default_controls.size();i<default_controls.size()+2;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==default_controls.size())
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.4*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
                }
                else if (i==default_controls.size()+1)
                {
                    _buttons[i]._text.setString("Default");
                    _buttons[i]._target="Default";
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.6*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));

                if (i==default_controls.size())
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.4*_sfac*raw_width),int(height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio)));
                }
                else if (i==default_controls.size()+1)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.6*_sfac*raw_width),int(height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio)));
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
        }
        void update (double dt,sf::Vector2i mouse_pos)
        {
            _shouldUpdate=false;

            //check for key conflicts.
            std::map<std::string,bool> result=getControlConflicts();

            sf::Color red1=sf::Color(170,0,0,150);
            sf::Color red2=sf::Color(170,0,0,150);
            sf::Color red3=sf::Color(170,0,0,200);

            int colour1[4]={100,100,100,150};
            int colour2[4]={100,100,100,150};
            int colour3[4]={169,169,169,200};

            for (auto it: result)
            {
                //make the button red.
                for (int i=0;i<_buttons.size();i++)
                {
                    if (_buttons[i]._target==it.first)
                    {
                        if (it.second==true)
                        {
                            _buttons[i]._colour1=red1;
                            _buttons[i]._colour2=red2;
                            _buttons[i]._colour3=red3;
                        }
                        else
                        {
                            _buttons[i]._colour1=sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]);
                            _buttons[i]._colour2=sf::Color(colour2[0],colour2[1],colour2[2],colour2[3]);
                            _buttons[i]._colour3=sf::Color(colour3[0],colour3[1],colour3[2],colour3[3]);
                        }
                    }
                }
            }
        };
};

#endif // CONTROL-SCREEN_H_INCLUDED
