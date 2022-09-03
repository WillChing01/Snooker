#ifndef CONTROL-SCREEN_H_INCLUDED
#define CONTROL-SCREEN_H_INCLUDED

class ControlScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;
        sf::Text universalControlsText;
        sf::Text shotControlsText;
        sf::Text placeControlsText;
        sf::Text infoText;

        sf::RectangleShape universalControlsShape;
        sf::RectangleShape shotControlsShape;
        sf::RectangleShape placeControlsShape;

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

            const int keyCharSize=int(_sfac*raw_height*0.02);
            const int controlShapesOffset=_sfac*raw_width*0.07;

            universalControlsText.setFont(_boldfont);
            universalControlsText.setString("Universal controls:");
            universalControlsText.setCharacterSize(keyCharSize);
            textrect=universalControlsText.getLocalBounds();
            universalControlsText.setOrigin(sf::Vector2f(0.,int(textrect.top+0.5*textrect.height)));
            universalControlsText.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5-controlShapesOffset),int(_sfac*raw_height*0.225)));
            universalControlsText.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&universalControlsText);

            shotControlsText.setFont(_boldfont);
            shotControlsText.setString("Shot controls:");
            shotControlsText.setCharacterSize(keyCharSize);
            textrect=shotControlsText.getLocalBounds();
            shotControlsText.setOrigin(sf::Vector2f(0.,int(textrect.top+0.5*textrect.height)));
            shotControlsText.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5-controlShapesOffset),int(_sfac*raw_height*0.25)));
            shotControlsText.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&shotControlsText);

            placeControlsText.setFont(_boldfont);
            placeControlsText.setString("Place ball controls:");
            placeControlsText.setCharacterSize(keyCharSize);
            textrect=placeControlsText.getLocalBounds();
            placeControlsText.setOrigin(sf::Vector2f(0.,int(textrect.top+0.5*textrect.height)));
            placeControlsText.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5-controlShapesOffset),int(_sfac*raw_height*0.275)));
            placeControlsText.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&placeControlsText);

            infoText.setFont(_thinfont);
            infoText.setString("Controls in the same category can overlap but universal controls can overlap with all others");
            infoText.setCharacterSize(int(keyCharSize*0.75));
            textrect=infoText.getLocalBounds();
            infoText.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            infoText.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.3)));
            infoText.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&infoText);

            universalControlsShape.setSize(sf::Vector2f(keyCharSize,keyCharSize));
            universalControlsShape.setFillColor(universalControlsColour);
            textrect=universalControlsShape.getLocalBounds();
            universalControlsShape.setOrigin(sf::Vector2f(textrect.left+textrect.width,textrect.top+0.5*textrect.height));
            textrect=universalControlsText.getGlobalBounds();
            universalControlsShape.setPosition(sf::Vector2f(_sfac*raw_width*0.5+controlShapesOffset,textrect.top+textrect.height*0.5));
            _shapes.push_back(&universalControlsShape);

            shotControlsShape.setSize(sf::Vector2f(keyCharSize,keyCharSize));
            shotControlsShape.setFillColor(shotControlsColour);
            textrect=shotControlsShape.getLocalBounds();
            shotControlsShape.setOrigin(sf::Vector2f(textrect.left+textrect.width,textrect.top+0.5*textrect.height));
            textrect=shotControlsText.getGlobalBounds();
            shotControlsShape.setPosition(sf::Vector2f(_sfac*raw_width*0.5+controlShapesOffset,textrect.top+textrect.height*0.5));
            _shapes.push_back(&shotControlsShape);

            placeControlsShape.setSize(sf::Vector2f(keyCharSize,keyCharSize));
            placeControlsShape.setFillColor(placeControlsColour);
            textrect=placeControlsShape.getLocalBounds();
            placeControlsShape.setOrigin(sf::Vector2f(textrect.left+textrect.width,textrect.top+0.5*textrect.height));
            textrect=placeControlsText.getGlobalBounds();
            placeControlsShape.setPosition(sf::Vector2f(_sfac*raw_width*0.5+controlShapesOffset,textrect.top+textrect.height*0.5));
            _shapes.push_back(&placeControlsShape);

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

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }
        }
        void update (double dt,sf::Vector2i mouse_pos)
        {
            _shouldUpdate=false;

            for (int i=0;i<_buttons.size();i++)
            {
                if (!_buttons[i]._isactive) {continue;}
                if (!_buttons[i]._shouldExecute) {continue;}

                _buttons[i]._shouldExecute=false;

                if (_buttons[i]._target=="Default")
                {
                    user_controls=default_controls;
                    sf::FloatRect bounds;
                    for (int j=0;j<default_controls.size();j++)
                    {
                        _buttons[j]._text.setString(KeyToString(user_controls[_buttons[j]._target]));
                        bounds=_buttons[j]._text.getLocalBounds();
                        _buttons[j]._text.setOrigin(sf::Vector2f(int(bounds.left+0.5*bounds.width),int(bounds.top+0.5*bounds.height)));
                    }

                    //update the config file with new controls.
                    std::ofstream file(_userConfigFile,std::ofstream::out | std::ofstream::trunc);
                    for (auto thing:user_controls) {file << KeyToString(thing.second) << "\n";}
                    file.close();
                }

                break;
            }

            //check for key conflicts.
            std::map<std::string,bool> result=getControlConflicts();

            sf::Color red1=sf::Color(170,0,0,150);
            sf::Color red2=sf::Color(170,0,0,150);
            sf::Color red3=sf::Color(170,0,0,200);

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
                            _buttons[i]._colour1=buttonColour1;
                            _buttons[i]._colour2=buttonColour2;
                            _buttons[i]._colour3=buttonColour3;
                        }
                    }
                }
            }
        };
};

#endif // CONTROL-SCREEN_H_INCLUDED
