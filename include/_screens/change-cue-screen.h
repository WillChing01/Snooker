#ifndef CHANGE-CUE-SCREEN_H_INCLUDED
#define CHANGE-CUE-SCREEN_H_INCLUDED

//button callbacks
void scrollLeft(GameState* game_state, std::map<std::string,std::string>* payload);
void scrollRight(GameState* game_state, std::map<std::string,std::string>* payload);
void selectCue(GameState* game_state, std::map<std::string,std::string>* payload);

class ChangeCueScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;
        sf::Text pageNumber;

        std::vector<sf::Sprite> cuepng;
        std::vector<sf::Texture> cuetexture;

        sf::RectangleShape selectionBorder;

        const double buttonwidth=_sfac*raw_width*0.15;

    public:
        int _page=0;
        const int totalCues=9;
        const int numPerPage=5;

        int selectedCueIndex=std::stoi(cuetexturefile.substr(_cueFilePrefix.length()+3).substr(0,cuetexturefile.length()-_cueFilePrefix.length()-3-4));

        ChangeCueScreen(double sfac=dfactor) : GameState(sfac)
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

            std::string custom_input;
            std::ifstream cuefile(_userCueConfigFile);
            if (cuefile.is_open())
            {
                while (getline(cuefile,custom_input))
                {
                    cuetexturefile=custom_input;
                } cuefile.close();
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Change cue");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            title.setPosition(sf::Vector2f(int(_sfac*raw_width*0.5),int(_sfac*raw_height*0.15)));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

            sf::FloatRect bounds;

            _buttons.push_back(RectButton());

            for (int i=0;i<numPerPage;i++)
            {
                cuepng.push_back(sf::Sprite());
                _buttons.push_back(RectButton());
            }

            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());

            for (int i=0;i<totalCues;i++)
            {
                cuetexture.push_back(sf::Texture());
                cuetexture[i].loadFromFile(_cueFilePrefix+"cue"+std::to_string(i)+".png");
                cuetexture[i].setSmooth(true);
            }

            for (int i=0;i<numPerPage;i++)
            {
                cuepng[i].setTexture(cuetexture[i]);
                cuepng[i].scale(65.*dfactor/5213.,65.*dfactor/5213.);
                bounds=cuepng[i].getLocalBounds();
                cuepng[i].setOrigin(bounds.left,bounds.top+bounds.height*0.5);
                cuepng[i].setPosition(sf::Vector2f(0.20*_sfac*raw_width,0.30*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[0]._ratio));
                _shapes.push_back(&cuepng[i]);
            }

            //set up the selection border.
            const double factor=1.05;
            bounds=cuepng[0].getLocalBounds();
            selectionBorder.setSize(sf::Vector2f(factor*bounds.width*65.*dfactor/5213.,buttonwidth/_buttons[0]._ratio));

            bounds=selectionBorder.getLocalBounds();
            selectionBorder.setOrigin(bounds.left+bounds.width*(1.-1./factor)*0.5,bounds.top+bounds.height*0.5);

            if (selectedCueIndex<numPerPage)
            {
                sf::Vector2f pos=cuepng[selectedCueIndex].getPosition();
                selectionBorder.setPosition(pos);
            }
            else
            {
                selectionBorder.setPosition(sf::Vector2f(-1000.,-1000.));
            }

            selectionBorder.setOutlineColor(buttonOutlineColour2);
            selectionBorder.setFillColor(sf::Color(0,0,0,0));
            selectionBorder.setOutlineThickness(absOutlineThickness);

            _shapes.push_back(&selectionBorder);

            for (int i=0;i<numPerPage+3;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (i==0)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==numPerPage+1 || i==numPerPage+2)
                {
                    _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth/_buttons[i]._ratio,buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._shape.setOrigin(0.5*buttonwidth/_buttons[i]._ratio,0.5*buttonwidth/_buttons[i]._ratio);

                    if (i==numPerPage+1)
                    {
                        _buttons[i]._shape.setPosition(sf::Vector2f(0.45*_sfac*raw_width,0.31*_sfac*raw_height+(numPerPage*1.2)*buttonwidth/_buttons[i]._ratio));
                    }
                    else
                    {
                        _buttons[i]._shape.setPosition(sf::Vector2f(0.55*_sfac*raw_width,0.31*_sfac*raw_height+(numPerPage*1.2)*buttonwidth/_buttons[i]._ratio));
                    }
                }
                else
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.75*_sfac*raw_width,0.30*_sfac*raw_height+((i-1)*1.2)*buttonwidth/_buttons[i]._ratio));
                }

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                }
                else if (i==numPerPage+1 || i==numPerPage+2)
                {
                    if (i==numPerPage+1)
                    {
                        _buttons[i]._text.setString("<");
                        _buttons[i]._callback=std::bind(scrollLeft,std::placeholders::_1,std::placeholders::_2);
                    }
                    else
                    {
                        _buttons[i]._text.setString(">");
                        _buttons[i]._callback=std::bind(scrollRight,std::placeholders::_1,std::placeholders::_2);
                    }
                }
                else
                {
                    _buttons[i]._text.setString("Select");
                    _buttons[i]._target="Select"+std::to_string(i-1);
                    _buttons[i]._payload["Select"]=std::to_string(i-1);
                    _buttons[i]._callback=std::bind(selectCue,std::placeholders::_1,std::placeholders::_2);
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));

                if (i==0)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.5*_sfac*raw_width),int(0.50*_sfac*raw_height+(4*1.2)*buttonwidth/_buttons[i]._ratio)));
                }
                else if (i==numPerPage+1 || i==numPerPage+2)
                {
                    if (i==numPerPage+1)
                    {
                        _buttons[i]._text.setPosition(sf::Vector2f(0.45*_sfac*raw_width,0.31*_sfac*raw_height+(numPerPage*1.2)*buttonwidth/_buttons[i]._ratio));
                    }
                    else
                    {
                        _buttons[i]._text.setPosition(sf::Vector2f(0.55*_sfac*raw_width,0.31*_sfac*raw_height+(numPerPage*1.2)*buttonwidth/_buttons[i]._ratio));
                    }
                }
                else
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(int(0.75*_sfac*raw_width),int(0.30*_sfac*raw_height+((i-1)*1.2)*buttonwidth/_buttons[i]._ratio)));
                }

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }

            pageNumber.setFont(_thinfont);
            pageNumber.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
            pageNumber.setString("1");
            textrect=pageNumber.getLocalBounds();
            pageNumber.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));
            pageNumber.setPosition(sf::Vector2f(0.50*_sfac*raw_width,0.31*_sfac*raw_height+(numPerPage*1.2)*buttonwidth/_buttons[0]._ratio));
            pageNumber.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&pageNumber);
        }

        void update(double dt,sf::Vector2i mouse_pos)
        {
            _shouldUpdate=false;

            pageNumber.setString(std::to_string(_page+1));
            sf::FloatRect textrect=pageNumber.getLocalBounds();
            pageNumber.setOrigin(sf::Vector2f(int(textrect.left+0.5*textrect.width),int(textrect.top+0.5*textrect.height)));

            for (int i=_page*numPerPage;i<_page*numPerPage+numPerPage;i++)
            {
                if (i>=totalCues)
                {
                    // not enough cues to fill the page, must move away the texture/button.
                    cuepng[i%numPerPage].setPosition(sf::Vector2f(-1000.,-1000.));
                    _buttons[i%numPerPage+1]._shape.setPosition(sf::Vector2f(-1000.,-1000.));
                    _buttons[i%numPerPage+1]._text.setPosition(sf::Vector2f(-1000.,-1000.));
                }
                else
                {
                    //update cue image.
                    cuepng[i%numPerPage].setTexture(cuetexture[i]);
                    sf::FloatRect bounds=cuepng[i%numPerPage].getLocalBounds();
                    cuepng[i%numPerPage].setOrigin(bounds.left,bounds.top+bounds.height*0.5);
                    cuepng[i%numPerPage].setPosition(sf::Vector2f(0.20*_sfac*raw_width,0.30*_sfac*raw_height+((i%numPerPage)*1.2)*buttonwidth/_buttons[i%numPerPage]._ratio));

                    //update button.
                    _buttons[i%numPerPage+1]._shape.setPosition(sf::Vector2f(0.75*_sfac*raw_width,0.30*_sfac*raw_height+((i%numPerPage)*1.2)*buttonwidth/_buttons[i%numPerPage]._ratio));
                    _buttons[i%numPerPage+1]._text.setPosition(sf::Vector2f(int(0.75*_sfac*raw_width),int(0.30*_sfac*raw_height+((i%numPerPage)*1.2)*buttonwidth/_buttons[i%numPerPage]._ratio)));
                    _buttons[i%numPerPage+1]._target="Select"+std::to_string(i);
                    _buttons[i%numPerPage+1]._payload["Select"]=std::to_string(i);
                }
            }

            if (selectedCueIndex>=_page*numPerPage && selectedCueIndex<_page*numPerPage+numPerPage)
            {
                sf::Vector2f pos=cuepng[selectedCueIndex%numPerPage].getPosition();
                selectionBorder.setPosition(pos);
            }
            else
            {
                selectionBorder.setPosition(sf::Vector2f(-1000.,-1000.));
            }
        };
};

#endif // CHANGE-CUE-SCREEN_H_INCLUDED
