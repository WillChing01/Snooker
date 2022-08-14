#ifndef CHANGE-CUE-SCREEN_H_INCLUDED
#define CHANGE-CUE-SCREEN_H_INCLUDED

class ChangeCueScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;
        sf::Sprite cuepng[3];
        sf::Texture cuetexture[3];

        sf::RectangleShape selectionBorder;

    public:
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

            int textcolour1[4]={255,255,255,255};
            int textcolour2[4]={255,255,255,255};
            int textcolour3[4]={255,255,255,255};
            int outlinecolour1[4]={200,200,200,255};
            int outlinecolour2[4]={255,255,255,255};
            int outlinecolour3[4]={255,255,255,255};
            int colour1[4]={100,100,100,150};
            int colour2[4]={100,100,100,150};
            int colour3[4]={169,169,169,200};

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Change cue");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

            sf::FloatRect bounds;
            double buttonwidth=_sfac*raw_width*0.15;

            for (int i=0;i<4;i++)
            {
                _buttons.push_back(RectButton());
            }

            for (int i=0;i<3;i++)
            {
                cuetexture[i].loadFromFile(_cueFilePrefix+"cue"+std::to_string(i)+".png");
                cuetexture[i].setSmooth(true);
                cuepng[i].setTexture(cuetexture[i]);
                cuepng[i].scale(65.*dfactor/5213.,65.*dfactor/5213.);
                bounds=cuepng[i].getLocalBounds();
                cuepng[i].setOrigin(bounds.left,bounds.top+bounds.height*0.5);
                cuepng[i].setPosition(sf::Vector2f(0.20*_sfac*raw_width,0.30*_sfac*raw_height+((i+1)*1.2)*buttonwidth/_buttons[0]._ratio));
                _shapes.push_back(&cuepng[i]);
            }

            //set up the selection border.
            const double factor=1.05;
            bounds=cuepng[0].getLocalBounds();
            selectionBorder.setSize(sf::Vector2f(factor*bounds.width*65.*dfactor/5213.,buttonwidth/_buttons[0]._ratio));

            int selectedCueIndex=cuetexturefile[cuetexturefile.size()-5]-'0';

            bounds=selectionBorder.getLocalBounds();
            selectionBorder.setOrigin(bounds.left+bounds.width*(1.-1./factor)*0.5,bounds.top+bounds.height*0.5);
            sf::Vector2f pos=cuepng[selectedCueIndex].getPosition();
            selectionBorder.setPosition(pos);

            selectionBorder.setOutlineColor(sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]));
            selectionBorder.setFillColor(sf::Color(0,0,0,0));
            selectionBorder.setOutlineThickness(absOutlineThickness);

            _shapes.push_back(&selectionBorder);

            for (int i=0;i<4;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (i==0)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.75*_sfac*raw_width,0.30*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                }
                else
                {
                    _buttons[i]._text.setString("Select");
                    _buttons[i]._target="Select"+std::to_string(i-1);
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));

                if (i==0)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.75*_sfac*raw_width,0.30*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
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

        void update(double dt,sf::Vector2i mouse_pos)
        {
            int selectedCueIndex=cuetexturefile[cuetexturefile.size()-5]-'0';

            sf::Vector2f pos=cuepng[selectedCueIndex].getPosition();
            selectionBorder.setPosition(pos);
        };
};

#endif // CHANGE-CUE-SCREEN_H_INCLUDED
