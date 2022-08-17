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

                if (i==0)
                {
                    t[i].setString("Aim left");
                    _buttons[i]._target="Aim left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Aim left"]));
                }
                else if (i==1)
                {
                    t[i].setString("Aim right");
                    _buttons[i]._target="Aim right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Aim right"]));
                }
                else if (i==2)
                {
                    t[i].setString("Precise aim left");
                    _buttons[i]._target="Precise aim left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Precise aim left"]));
                }
                else if (i==3)
                {
                    t[i].setString("Precise aim right");
                    _buttons[i]._target="Precise aim right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Precise aim right"]));
                }
                else if (i==4)
                {
                    t[i].setString("Increase power");
                    _buttons[i]._target="Increase power";
                    _buttons[i]._text.setString(KeyToString(user_controls["Increase power"]));
                }
                else if (i==5)
                {
                    t[i].setString("Decrease power");
                    _buttons[i]._target="Decrease power";
                    _buttons[i]._text.setString(KeyToString(user_controls["Decrease power"]));
                }
                else if (i==6)
                {
                    t[i].setString("Increase cue elevation");
                    _buttons[i]._target="Increase cue elevation";
                    _buttons[i]._text.setString(KeyToString(user_controls["Increase cue elevation"]));
                }
                else if (i==7)
                {
                    t[i].setString("Decrease cue elevation");
                    _buttons[i]._target="Decrease cue elevation";
                    _buttons[i]._text.setString(KeyToString(user_controls["Decrease cue elevation"]));
                }
                else if (i==8)
                {
                    t[i].setString("Offset cue tip up");
                    _buttons[i]._target="Offset cue tip up";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip up"]));
                }
                else if (i==9)
                {
                    t[i].setString("Offset cue tip down");
                    _buttons[i]._target="Offset cue tip down";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip down"]));
                }
                else if (i==10)
                {
                    t[i].setString("Offset cue tip left");
                    _buttons[i]._target="Offset cue tip left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip left"]));
                }
                else if (i==11)
                {
                    t[i].setString("Offset cue tip right");
                    _buttons[i]._target="Offset cue tip right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip right"]));
                }
                else if (i==12)
                {
                    t[i].setString("Strike cueball");
                    _buttons[i]._target="Strike cueball";
                    _buttons[i]._text.setString(KeyToString(user_controls["Strike cueball"]));
                }
                else if (i==13)
                {
                    t[i].setString("Move ball up");
                    _buttons[i]._target="Move ball up";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball up"]));
                }
                else if (i==14)
                {
                    t[i].setString("Move ball down");
                    _buttons[i]._target="Move ball down";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball down"]));
                }
                else if (i==15)
                {
                    t[i].setString("Move ball left");
                    _buttons[i]._target="Move ball left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball left"]));
                }
                else if (i==16)
                {
                    t[i].setString("Move ball right");
                    _buttons[i]._target="Move ball right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball right"]));
                }
                else if (i==17)
                {
                    t[i].setString("Place ball");
                    _buttons[i]._target="Place ball";
                    _buttons[i]._text.setString(KeyToString(user_controls["Place ball"]));
                }
                else if (i==18)
                {
                    t[i].setString("Toggle mute");
                    _buttons[i]._target="Toggle mute";
                    _buttons[i]._text.setString(KeyToString(user_controls["Toggle mute"]));
                }
                else if (i==19)
                {
                    t[i].setString("Pause game");
                    _buttons[i]._target="Pause game";
                    _buttons[i]._text.setString(KeyToString(user_controls["Pause game"]));
                }

                textrect=t[i].getLocalBounds();
                t[i].setOrigin(sf::Vector2f(int(textrect.left),int(textrect.top)));
                t[i].setFillColor(sf::Color(255,255,255));


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

            sf::Color red1=sf::Color(170,0,0);
            sf::Color red2=sf::Color(170,0,0);
            sf::Color red3=sf::Color(200,0,0);

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
