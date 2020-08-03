#ifndef CLIENT_H_INCLUDED
#define CLIENT_H_INCLUDED

#include <string>
#include <map>
#include <SFML/Graphics.hpp>
#include <SFML/Network.hpp>
#include "objects.h"

const std::map<std::string,sf::Keyboard::Key> default_controls=
{
{"Aim left",sf::Keyboard::Left},
{"Aim right",sf::Keyboard::Right},
{"Precise aim left",sf::Keyboard::Comma},
{"Precise aim right",sf::Keyboard::Period},
{"Increase power",sf::Keyboard::Up},
{"Decrease power",sf::Keyboard::Down},
{"Increase cue elevation",sf::Keyboard::LBracket},
{"Decrease cue elevation",sf::Keyboard::RBracket},
{"Offset cue tip up",sf::Keyboard::Numpad8},
{"Offset cue tip down",sf::Keyboard::Numpad2},
{"Offset cue tip left",sf::Keyboard::Numpad4},
{"Offset cue tip right",sf::Keyboard::Numpad6},
{"Strike cueball",sf::Keyboard::Space},
{"Move ball up",sf::Keyboard::Up},
{"Move ball down",sf::Keyboard::Down},
{"Move ball left",sf::Keyboard::Left},
{"Move ball right",sf::Keyboard::Right},
{"Place ball",sf::Keyboard::Enter},
{"Pause game",sf::Keyboard::Escape}
};

std::map<std::string,sf::Keyboard::Key> user_controls=default_controls;

#define ITEM(x) case sf::Keyboard::x : return #x;

std::string KeyToString(sf::Keyboard::Key k)
{
    switch(k)
    {
        ITEM(A);
        ITEM(B);
        ITEM(C);
        ITEM(D);
        ITEM(E);
        ITEM(F);
        ITEM(G);
        ITEM(H);
        ITEM(I);
        ITEM(J);
        ITEM(K);
        ITEM(L);
        ITEM(M);
        ITEM(N);
        ITEM(O);
        ITEM(P);
        ITEM(Q);
        ITEM(R);
        ITEM(S);
        ITEM(T);
        ITEM(U);
        ITEM(V);
        ITEM(W);
        ITEM(X);
        ITEM(Y);
        ITEM(Z);
        ITEM(Num0);
        ITEM(Num1);
        ITEM(Num2);
        ITEM(Num3);
        ITEM(Num4);
        ITEM(Num5);
        ITEM(Num6);
        ITEM(Num7);
        ITEM(Num8);
        ITEM(Num9);
        ITEM(Escape);
        ITEM(LControl);
        ITEM(LShift);
        ITEM(LAlt);
        ITEM(LSystem);
        ITEM(RControl);
        ITEM(RShift);
        ITEM(RAlt);
        ITEM(RSystem);
        ITEM(Menu);
        ITEM(LBracket);
        ITEM(RBracket);
        ITEM(Semicolon);
        ITEM(Comma);
        ITEM(Period);
        ITEM(Quote);
        ITEM(Slash);
        ITEM(Backslash);
        ITEM(Tilde);
        ITEM(Equal);
        ITEM(Hyphen);
        ITEM(Space);
        ITEM(Enter);
        ITEM(Backspace);
        ITEM(Tab);
        ITEM(PageUp);
        ITEM(PageDown);
        ITEM(End);
        ITEM(Home);
        ITEM(Insert);
        ITEM(Delete);
        ITEM(Add);
        ITEM(Subtract);
        ITEM(Multiply);
        ITEM(Divide);
        ITEM(Left);
        ITEM(Right);
        ITEM(Up);
        ITEM(Down);
        ITEM(Numpad0);
        ITEM(Numpad1);
        ITEM(Numpad2);
        ITEM(Numpad3);
        ITEM(Numpad4);
        ITEM(Numpad5);
        ITEM(Numpad6);
        ITEM(Numpad7);
        ITEM(Numpad8);
        ITEM(Numpad9);
        ITEM(F1);
        ITEM(F2);
        ITEM(F3);
        ITEM(F4);
        ITEM(F5);
        ITEM(F6);
        ITEM(F7);
        ITEM(F8);
        ITEM(F9);
        ITEM(F10);
        ITEM(F11);
        ITEM(F12);
        ITEM(F13);
        ITEM(F14);
        ITEM(F15);
        ITEM(Pause);
        ITEM(KeyCount);
        ITEM(Unknown);
    }
}

class RectButton
{
    public:
        double _ratio=4.;
        //ratio is width/height.
        double _textfactor=0.75;
        double _absoutlinethickness=2.;

        bool _controlchange=false;

        std::string _target;

        sf::RectangleShape _shape;
        sf::Font _font;
        sf::Text _text;

        sf::Color _textcolour1;
        sf::Color _textcolour2;
        sf::Color _textcolour3;

        sf::Color _colour1;
        sf::Color _colour2;
        sf::Color _colour3;

        sf::Color _outlinecolour1;
        sf::Color _outlinecolour2;
        sf::Color _outlinecolour3;

        RectButton() {};
};

class InputBox
{
    public:
        double _textfactor=0.75;
        double _absoutlinethickness=2.;
        double _abscursorthickness=1.;

        double _t=0.;
        double _ton=0.8;
        double _toff=0.8;

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

class GameState
{
    public:
        std::vector<RectButton> _buttons;
        std::vector<InputBox> _inputboxes;
        std::vector<sf::Drawable*> _shapes;

        sf::Color _background;

        double _sfac;

        GameState(double sfac) {_sfac=sfac;}
        virtual void update(double dt)=0;
};

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

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
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

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
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
        void update(double dt) {};
};

class OptionsScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;

    public:
        OptionsScreen(double sfac=dfactor) : GameState(sfac)
        {
            _background=sf::Color(0,0,0);

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Options");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
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

            for (int i=0;i<3;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Controls");
                    _buttons[i]._target="Controls";
                }
                else if (i==1)
                {
                    _buttons[i]._text.setString("Change cue");
                    _buttons[i]._target="Changecue";
                }
                else if (i==2)
                {
                    _buttons[i]._text.setString("Back");
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
        void update(double dt) {};
};

class ControlScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;

        sf::Text t[18];

    public:
        ControlScreen(double sfac=dfactor) : GameState(sfac)
        {
            _background=sf::Color(0,0,0);

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Controls");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());

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

            for (int i=0;i<18;i++)
            {
                t[i].setFont(_thinfont);
                t[i].setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));

                _buttons[i]._controlchange=true;
                _buttons[i]._shape.setSize(sf::Vector2f(4.*buttonwidth/_buttons[i]._ratio,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(2.*buttonwidth/_buttons[i]._ratio,0.5*buttonwidth/_buttons[i]._ratio);

                if (i<9)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.475*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.775*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+((i-9)*1.2)*buttonwidth/_buttons[i]._ratio));
                }

                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
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
                else if (i==4)
                {
                    t[i].setString("Decrease power");
                    _buttons[i]._target="Decrease power";
                    _buttons[i]._text.setString(KeyToString(user_controls["Decrease power"]));
                }
                else if (i==5)
                {
                    t[i].setString("Increase cue elevation");
                    _buttons[i]._target="Increase cue elevation";
                    _buttons[i]._text.setString(KeyToString(user_controls["Increase cue elevation"]));
                }
                else if (i==6)
                {
                    t[i].setString("Decrease cue elevation");
                    _buttons[i]._target="Decrease cue elevation";
                    _buttons[i]._text.setString(KeyToString(user_controls["Decrease cue elevation"]));
                }
                else if (i==7)
                {
                    t[i].setString("Offset cue tip up");
                    _buttons[i]._target="Offset cue tip up";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip up"]));
                }
                else if (i==8)
                {
                    t[i].setString("Offset cue tip down");
                    _buttons[i]._target="Offset cue tip down";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip down"]));
                }
                else if (i==9)
                {
                    t[i].setString("Offset cue tip left");
                    _buttons[i]._target="Offset cue tip left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip left"]));
                }
                else if (i==10)
                {
                    t[i].setString("Offset cue tip right");
                    _buttons[i]._target="Offset cue tip right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Offset cue tip right"]));
                }
                else if (i==11)
                {
                    t[i].setString("Strike cueball");
                    _buttons[i]._target="Strike cueball";
                    _buttons[i]._text.setString(KeyToString(user_controls["Strike cueball"]));
                }
                else if (i==12)
                {
                    t[i].setString("Move ball up");
                    _buttons[i]._target="Move ball up";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball up"]));
                }
                else if (i==13)
                {
                    t[i].setString("Move ball down");
                    _buttons[i]._target="Move ball down";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball down"]));
                }
                else if (i==14)
                {
                    t[i].setString("Move ball left");
                    _buttons[i]._target="Move ball left";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball left"]));
                }
                else if (i==15)
                {
                    t[i].setString("Move ball right");
                    _buttons[i]._target="Move ball right";
                    _buttons[i]._text.setString(KeyToString(user_controls["Move ball right"]));
                }
                else if (i==16)
                {
                    t[i].setString("Place ball");
                    _buttons[i]._target="Place ball";
                    _buttons[i]._text.setString(KeyToString(user_controls["Place ball"]));
                }
                else if (i==17)
                {
                    t[i].setString("Pause game");
                    _buttons[i]._target="Pause game";
                    _buttons[i]._text.setString(KeyToString(user_controls["Pause game"]));
                }

                textrect=t[i].getLocalBounds();
                t[i].setOrigin(sf::Vector2f(textrect.left,textrect.top));
                t[i].setFillColor(sf::Color(255,255,255));


                textrect=_buttons[i]._text.getLocalBounds();

                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));

                if (i<9)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.475*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                    t[i].setPosition(sf::Vector2f(0.225*_sfac*raw_width,height+(i*1.2)*buttonwidth/_buttons[i]._ratio-0.5*buttonwidth/_buttons[i]._ratio+0.5*0.25*buttonwidth/_buttons[i]._ratio));
                }
                else
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.775*_sfac*raw_width-2.*buttonwidth/_buttons[i]._ratio,height+((i-9)*1.2)*buttonwidth/_buttons[i]._ratio));
                    t[i].setPosition(sf::Vector2f(0.525*_sfac*raw_width,height+((i-9)*1.2)*buttonwidth/_buttons[i]._ratio-0.5*buttonwidth/_buttons[i]._ratio+0.5*0.25*buttonwidth/_buttons[i]._ratio));
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
            for (int i=18;i<20;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==18)
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.4*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
                }
                else if (i==19)
                {
                    _buttons[i]._text.setString("Default");
                    _buttons[i]._target="Default";
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.6*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));

                if (i==18)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.4*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
                }
                else if (i==19)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.6*_sfac*raw_width,height+(11*1.24)*(_sfac*raw_width*0.083)/_buttons[i]._ratio));
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
        void update (double dt) {};
};

class ChangeCueScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;

    public:
        ChangeCueScreen(double sfac=dfactor) : GameState(sfac)
        {
            _background=sf::Color(0,0,0);

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Change cue");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
            title.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&title);

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

            for (int i=0;i<1;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Back");
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
        void update(double dt) {};
};

class SingleplayerScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;

        sf::Text title;

    public:
        SingleplayerScreen(double sfac=dfactor) : GameState(sfac)
        {
            _background=sf::Color(0,0,0);

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Singleplayer");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
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

            for (int i=0;i<3;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Play vs. AI");
                    _buttons[i]._target="SingleplayerAI";
                }
                else if (i==1)
                {
                    _buttons[i]._text.setString("Solo lineup");
                    _buttons[i]._target="SingleplayerLineup";
                }
                else if (i==2)
                {
                    _buttons[i]._text.setString("Back");
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
        void update(double dt) {};
};

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

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            title.setFont(_boldfont);
            title.setCharacterSize(int(_sfac*raw_height*0.1));
            title.setString("Multiplayer");
            textrect=title.getLocalBounds();
            title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15));
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
            ip.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            ip.setPosition(sf::Vector2f(_sfac*raw_width*0.3,_sfac*raw_height*0.5));
            ip.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&ip);

            port.setFont(_thinfont);
            port.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
            port.setString("Port : "+portnum);
            textrect=port.getLocalBounds();
            port.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            port.setPosition(sf::Vector2f(_sfac*raw_width*0.3,_sfac*raw_height*0.5+1.2*buttonwidth/_buttons[0]._ratio));
            port.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&port);

            for (int i=0;i<3;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Host game");
                    _buttons[i]._target="MultiplayerHost";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.3*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(0.3*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==1)
                {
                    _buttons[i]._text.setString("Join game");
                    _buttons[i]._target="MultiplayerJoin";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.7*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(0.7*_sfac*raw_width,0.50*_sfac*raw_height+(2.*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==2)
                {
                    _buttons[i]._text.setString("Back");
                    _buttons[i]._target="Quit";
                    textrect=_buttons[i]._text.getLocalBounds();
                    _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4.*1.2)*buttonwidth/_buttons[i]._ratio));
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.50*_sfac*raw_height+(4.*1.2)*buttonwidth/_buttons[i]._ratio));
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

                if (!_inputboxes[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _inputboxes[i]._text.setFont(_inputboxes[i]._font);
                _inputboxes[i]._text.setCharacterSize(int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio));
                _inputboxes[i]._text.setFillColor(sf::Color(255,255,255));
                _inputboxes[i]._text.setString("");
                _inputboxes[i]._text.setPosition(sf::Vector2f(0.7*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio-0.5*int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio)));

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
                _inputboxes[i]._backtext.setOrigin(sf::Vector2f(textrect.left,textrect.top+textrect.height*0.5));
                _inputboxes[i]._backtext.setFillColor(sf::Color(255,255,255,150));
                _inputboxes[i]._backtext.setPosition(sf::Vector2f(0.7*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5+i*1.2*buttonwidth/_buttons[0]._ratio));

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
                    _inputboxes[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio-0.5*int(buttonwidth*_buttons[0]._textfactor/_buttons[0]._ratio)));
                    _inputboxes[i]._backtext.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio));
                    _inputboxes[i]._cursor.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.5*boxwidth+2.*_inputboxes[i]._absoutlinethickness,_sfac*raw_height*0.5-i*1.2*buttonwidth/_buttons[0]._ratio));
                }
            }
        }
        void update(double dt)
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

class GameScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;
        int gametype=0;
        //gametype 0 - multiplayer host.
        //gametype 1 - multiplayer but not host.
        //gametype 2 - singleplayer vs. AI.
        //gametype 3 - singleplayer lineup.
        Server server;
        std::string targetip;
        unsigned short targetport;
        sf::TcpSocket socket;
        sf::Packet packet;
        sf::Uint16 packetId=0;
        sf::Uint32 resultsize=0;
        std::array<double,66> temp;

        sf::Font _font;
        sf::Font _scorefont;
        std::string p1name="PLAYER 1";
        std::string p2name="PLAYER 2";
        int p1score=0;
        int p2score=0;
        int p1frames=0;
        int p2frames=0;
        int framesbestof=35;
        bool placing_white=true;
        bool touching=false;
        bool isyourturn=true;
        bool isfoul=false;

        bool change=true;
        bool done=true;
        int t=0;
        double dt=0.;

        double power=50.;
        double dist;
        double angle;

        bool ispaused=false;
        double pausea;
        bool ispausepressed=false;
        sf::RectangleShape pauserect;
        sf::Text pausetext;

        bool gameover=false;
        sf::RectangleShape gameoverrect;
        sf::Text gameovertext;

        //stats.
        sf::Text stats_title;
        std::vector<sf::Text> stats_text;
        std::vector<double> stats_text_y;
        int p1_highbreak=0;
        int p2_highbreak=0;
        int p1_centuries=0;
        int p2_centuries=0;


        sf::Vector2f spin_selector_pos;
        sf::Vector2f spin_dot_pos;

        sf::RectangleShape panel_background;
        sf::CircleShape spin_selector;
        sf::CircleShape spin_dot;
        sf::RectangleShape power_bar[100]={};
        sf::RectangleShape power_outline;
        sf::CircleShape elevation_border;
        sf::CircleShape elevation_ball;
        sf::RectangleShape elevation_pointer;
        sf::Text elevation_display;
        sf::Text textframesbestof;
        sf::Text textp1frames;
        sf::Text textp2frames;
        sf::Text textp1score;
        sf::Text textp2score;
        sf::Text textp1name;
        sf::Text textp2name;
        sf::CircleShape p1pointer;
        sf::CircleShape p2pointer;
        sf::RectangleShape framescoresrect;
        sf::RectangleShape p1scorerect;
        sf::RectangleShape p2scorerect;
        sf::RectangleShape p1namerect;
        sf::RectangleShape p2namerect;

        Cue cue;
        Ball balls[22];
        Cushion cushions[6];
        sf::CircleShape spots[6];
        sf::VertexArray baulkline;
        sf::VertexArray baulkcircle;

        Eigen::MatrixXd test=Eigen::MatrixXd::Constant(46,1,0.0);
        std::vector<std::array<double,66> > result;
        std::array<std::vector<double>,3> predict;

        sf::VertexArray cuetraj;
        sf::VertexArray cuetraj2;
        sf::VertexArray obtraj;
        sf::CircleShape ghostball;
        sf::CircleShape ghostball2;

    public:
        //set up the objects.
        GameScreen(double sfac=dfactor,int kind=0,std::string ip="", unsigned short port=50000,std::string name="PLAYER 1") : GameState(sfac)
        {
            gametype=kind;

            p1name=name;

            targetip=ip;
            targetport=port;

            //set up some stuff depending on the gametype here.

            _background=sf::Color(baizecolour[0],baizecolour[1],baizecolour[2]);

            if (!_thinfont.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            //set up pause rect.
            pauserect.setSize(sf::Vector2f(_sfac*raw_width,_sfac*raw_height));
            pauserect.setPosition(sf::Vector2f(0.,-_sfac*raw_height));
            pauserect.setFillColor(sf::Color(100,100,100,150));

            pausetext.setFont(_thinfont);
            pausetext.setCharacterSize(int(_sfac*raw_height*0.1));
            pausetext.setString("Paused");
            textrect=pausetext.getLocalBounds();
            pausetext.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            pausetext.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15-_sfac*raw_height));
            pausetext.setFillColor(sf::Color(255,255,255));

            //set up game over screen.
            gameoverrect.setSize(sf::Vector2f(_sfac*raw_width,_sfac*raw_height));
            gameoverrect.setPosition(sf::Vector2f(0.,-_sfac*raw_height));
            gameoverrect.setFillColor(sf::Color(100,100,100,150));

            gameovertext.setFont(_thinfont);
            gameovertext.setCharacterSize(int(_sfac*raw_height*0.1));
            gameovertext.setString("Game over");
            textrect=gameovertext.getLocalBounds();
            gameovertext.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            gameovertext.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15-_sfac*raw_height));
            gameovertext.setFillColor(sf::Color(255,255,255));

            int sh=int(_sfac*raw_height*0.04);

            stats_title.setFont(_thinfont);
            stats_title.setCharacterSize(int(sh*1.2));
            stats_title.setString("Overall Stats");
            textrect=stats_title.getLocalBounds();
            stats_title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            stats_title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.25-_sfac*raw_height));
            stats_title.setFillColor(sf::Color(255,255,255));

            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());

            double starth=0.35*_sfac*raw_height;
            double soffset=0.15*_sfac*raw_width;

            for (int i=0;i<stats_text.size();i++)
            {
                stats_text[i].setFont(_thinfont);
                stats_text[i].setCharacterSize(sh);

                if (i==0) {stats_text[i].setString("versus");}
                else if (i==1) {stats_text[i].setString(p1name);}
                else if (i==2) {stats_text[i].setString(p2name);}
                else if (i==3) {stats_text[i].setString("Frames");}
                else if (i==4) {stats_text[i].setString(std::to_string(p1frames));}
                else if (i==5) {stats_text[i].setString(std::to_string(p2frames));}
                else if (i==6) {stats_text[i].setString("High break");}
                else if (i==7) {stats_text[i].setString(std::to_string(p1_highbreak));}
                else if (i==8) {stats_text[i].setString(std::to_string(p2_highbreak));}
                else if (i==9) {stats_text[i].setString("Centuries");}
                else if (i==10) {stats_text[i].setString(std::to_string(p1_centuries));}
                else if (i==11) {stats_text[i].setString(std::to_string(p2_centuries));}

                textrect=stats_text[i].getLocalBounds();
                stats_text[i].setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                if (i%3==0)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                else if (i%3==1)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width-soffset,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                else if (i%3==2)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width+soffset,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                stats_text[i].setFillColor(sf::Color(255,255,255));

                stats_text_y.push_back(stats_text[i].getPosition().y);
            }
//            sf::Text stats_title;
//            sf::Text stats_p1name;
//            sf::Text stats_p2name;
//            sf::Text stats_highbreaktext;
//            sf::Text stats_centuriestext;
//            sf::Text stats_framestext;
//            sf::Text stats_p1frames;
//            sf::Text stats_p2frames;
//            sf::Text stats_p1highbreak;
//            sf::Text stats_p2highbreak;
//            sf::Text stats_p1centuries;
//            sf::Text stats_p2centuries;

            pausea=2.*_sfac*raw_height;

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

            for (int i=0;i<2;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                if (i==0)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==1)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.75*_sfac*raw_height));
                }

                if (!_buttons[i]._font.loadFromFile("Roboto-Thin.ttf")) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Exit");
                    _buttons[i]._target="Quit";
                }
                if (i==1)
                {
                    _buttons[i]._text.setString("Exit");
                    _buttons[i]._target="Quit";
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));

                if (i==0)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==1)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.75*_sfac*raw_height));
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

            //set up panel.
            panel_background=sf::RectangleShape(sf::Vector2f((_sfac*raw_width),(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            panel_background.setPosition(sf::Vector2f(0.,(_sfac*raw_height/(1.+panel_ratio))));
            panel_background.setFillColor(sf::Color(0,0,0));

            spin_selector.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            spin_selector.setOrigin((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            spin_selector.setPosition(sf::Vector2f(0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            spin_selector.setFillColor(sf::Color(255,255,255));
            spin_selector.setPointCount(200);

            spin_dot.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625);
            spin_dot.setOrigin(sf::Vector2f((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625));
            spin_dot.setPosition(sf::Vector2f(0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            spin_dot.setFillColor(sf::Color(255,0,0));
            spin_dot.setPointCount(100);

            for (int i=0;i<100;i++)
            {
                power_bar[i]=sf::RectangleShape(sf::Vector2f(0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.008*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
                power_bar[i].setPosition(sf::Vector2f(1.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+i*0.008*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
                power_bar[i].setFillColor(sf::Color(255,int(floor(2.55*i)),0));
            }
            sf::Color c;
            for (int i=0;i<100;i++)
            {
                c=power_bar[99-i].getFillColor();
                if (i<int(floor(power))) {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,255));}
                else {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,0));}
            }
            double outline=0.01*(_sfac*raw_height*panel_ratio/(1.+panel_ratio));
            power_outline=sf::RectangleShape(sf::Vector2f(0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+2*outline,0.8*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+2*outline));
            power_outline.setPosition(sf::Vector2f(1.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))-outline,(_sfac*raw_height/(1.+panel_ratio))+0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))-outline));
            power_outline.setFillColor(sf::Color(0,0,0));
            power_outline.setOutlineThickness(outline);
            power_outline.setOutlineColor(sf::Color(255,255,255));

            elevation_border.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            elevation_border.setOrigin((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            elevation_border.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_border.setFillColor(sf::Color(0,0,0));
            elevation_border.setOutlineThickness(outline);
            elevation_border.setOutlineColor(sf::Color(255,255,255));
            elevation_border.setPointCount(200);

            elevation_ball.setRadius(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)));
            elevation_ball.setOrigin(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)));
            elevation_ball.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_ball.setFillColor(sf::Color(255,255,255));
            elevation_ball.setPointCount(150);

            elevation_pointer=sf::RectangleShape(sf::Vector2f(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.005/0.02625));
            elevation_pointer.setOrigin(0,0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.005/0.02625);
            elevation_pointer.setPosition(sf::Vector2f((_sfac*raw_width)-0.3*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_pointer.setFillColor(sf::Color(255,0,0));

            if (!_font.loadFromFile("Roboto-Thin.ttf"))
            {
                std::cout << "ERROR. Couldn't load font!" << std::endl;
            }
            elevation_display.setFont(_font);
            elevation_display.setString("00°");
            elevation_display.setCharacterSize(int(0.15*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_display.setFillColor(sf::Color(0,0,0));
            sf::FloatRect textRect=elevation_display.getLocalBounds();
            elevation_display.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            elevation_display.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));

            double sbwidth=(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*6.;
            double sbheight=(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.2;

            if (!_scorefont.loadFromFile("Roboto-Bold.ttf"))
            {
                std::cout << "ERROR. Couldn't load font!" << std::endl;
            }
            textframesbestof.setFont(_scorefont);
            textframesbestof.setString("("+std::to_string(framesbestof)+")");
            textframesbestof.setCharacterSize(int(sbheight*0.7));
            textframesbestof.setFillColor(sf::Color(255,255,255));
            textRect=textframesbestof.getLocalBounds();
            textframesbestof.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textframesbestof.setPosition(sf::Vector2f((_sfac*raw_width)*0.5,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1frames.setFont(_scorefont);
            textp1frames.setString(std::to_string(p1frames));
            textp1frames.setCharacterSize(int(sbheight*0.7));
            textp1frames.setFillColor(sf::Color(255,255,255));
            textRect=textp1frames.getLocalBounds();
            textp1frames.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp1frames.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.06*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2frames.setFont(_scorefont);
            textp2frames.setString(std::to_string(p2frames));
            textp2frames.setCharacterSize(int(sbheight*0.7));
            textp2frames.setFillColor(sf::Color(255,255,255));
            textRect=textp2frames.getLocalBounds();
            textp2frames.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp2frames.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.06*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1score.setFont(_scorefont);
            textp1score.setString(std::to_string(p1score));
            textp1score.setCharacterSize(int(sbheight*0.7));
            textp1score.setFillColor(sf::Color(0,0,0));
            textRect=textp1score.getLocalBounds();
            textp1score.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp1score.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.1*sbwidth-0.5*0.15*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2score.setFont(_scorefont);
            textp2score.setString(std::to_string(p2score));
            textp2score.setCharacterSize(int(sbheight*0.7));
            textp2score.setFillColor(sf::Color(0,0,0));
            textRect=textp2score.getLocalBounds();
            textp2score.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp2score.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.1*sbwidth+0.5*0.15*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1name.setFont(_scorefont);
            textp1name.setString(p1name);
            textp1name.setCharacterSize(int(sbheight*0.7));
            textp1name.setFillColor(sf::Color(0,0,0));
            textRect=textp1name.getLocalBounds();
            textp1name.setOrigin(sf::Vector2f(0.,textRect.top+textRect.height/2.));
            textp1name.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.47*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2name.setFont(_scorefont);
            textp2name.setString(p2name);
            textp2name.setCharacterSize(int(sbheight*0.7));
            textp2name.setFillColor(sf::Color(0,0,0));
            textRect=textp2name.getLocalBounds();
            textp2name.setOrigin(sf::Vector2f(textRect.left+textRect.width,textRect.top+textRect.height/2.));
            textp2name.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.47*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            p1pointer.setRadius(sbheight*0.2);
            p1pointer.setPointCount(3);
            p1pointer.setFillColor(sf::Color(0,0,0));
            p1pointer.setOrigin(sbheight*0.2,sbheight*0.2);
            p1pointer.setRotation(-90.);
            p1pointer.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.105*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            p2pointer.setRadius(sbheight*0.2);
            p2pointer.setPointCount(3);
            p2pointer.setFillColor(sf::Color(0,0,0,0));
            p2pointer.setOrigin(sbheight*0.2,sbheight*0.2);
            p2pointer.setRotation(90.);
            p2pointer.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.105*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            framescoresrect=sf::RectangleShape(sf::Vector2f(sbwidth*0.2,sbheight));
            framescoresrect.setOrigin(sbwidth*0.1,0.5*sbheight);
            framescoresrect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            framescoresrect.setFillColor(sf::Color(51,153,255));

            p1scorerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.15,sbheight));
            p1scorerect.setOrigin(sbwidth*0.15*0.5,0.5*sbheight);
            p1scorerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.1*sbwidth-0.15*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p1scorerect.setFillColor(sf::Color(255,255,255));

            p2scorerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.15,sbheight));
            p2scorerect.setOrigin(sbwidth*0.15*0.5,0.5*sbheight);
            p2scorerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.1*sbwidth+0.15*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p2scorerect.setFillColor(sf::Color(255,255,255));

            p1namerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.25,sbheight));
            p1namerect.setOrigin(sbwidth*0.25*0.5,0.5*sbheight);
            p1namerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.25*sbwidth-0.25*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p1namerect.setFillColor(sf::Color(236,228,0));

            p2namerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.25,sbheight));
            p2namerect.setOrigin(sbwidth*0.25*0.5,0.5*sbheight);
            p2namerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.25*sbwidth+0.25*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p2namerect.setFillColor(sf::Color(236,228,0));

            _shapes.push_back(&panel_background);
            _shapes.push_back(&spin_selector);
            _shapes.push_back(&spin_dot);
            _shapes.push_back(&power_outline);
            for (int i=0;i<100;i++)
            {
                _shapes.push_back(&power_bar[i]);
            }
            _shapes.push_back(&elevation_border);
            _shapes.push_back(&elevation_ball);
            _shapes.push_back(&elevation_pointer);
            _shapes.push_back(&framescoresrect);
            _shapes.push_back(&p1scorerect);
            _shapes.push_back(&p2scorerect);
            _shapes.push_back(&p1namerect);
            _shapes.push_back(&p2namerect);
            _shapes.push_back(&p1pointer);
            _shapes.push_back(&p2pointer);
            _shapes.push_back(&elevation_display);
            _shapes.push_back(&textframesbestof);
            _shapes.push_back(&textp1frames);
            _shapes.push_back(&textp2frames);
            _shapes.push_back(&textp1score);
            _shapes.push_back(&textp2score);
            _shapes.push_back(&textp1name);
            _shapes.push_back(&textp2name);

            cushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
            cushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
            cushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
            cushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
            cushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
            cushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&cushions[i]._pocketshape);
            }
            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&cushions[i]._shape);
                _shapes.push_back(&cushions[i]._railshape);
                _shapes.push_back(&cushions[i]._p1shape);
                _shapes.push_back(&cushions[i]._p2shape);
            }

            for (int i=0;i<6;i++)
            {
                spots[i].setRadius(spot_r*_sfac);
                spots[i].setOrigin(spot_r*_sfac,spot_r*_sfac);
                spots[i].setFillColor(sf::Color(255,255,255));
            }

            for (int i=0;i<22;i++)
            {
                balls[i]._shape.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            }

            //cueball.
            balls[0]._shape.setFillColor(sf::Color(255,255,255));
            balls[0]._x=cueball_break_x;
            balls[0]._y=cueball_break_y;
            balls[0]._order=1;
            //yellow.
            balls[1]._shape.setFillColor(sf::Color(255,255,0));
            balls[1]._x=yellow_x;
            balls[1]._y=yellow_y;
            spots[0].setPosition(sf::Vector2f(yellow_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-yellow_y*_sfac));
            balls[1]._order=2;
            //green.
            balls[2]._shape.setFillColor(sf::Color(0,150,0));
            balls[2]._x=green_x;
            balls[2]._y=green_y;
            spots[1].setPosition(sf::Vector2f(green_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-green_y*_sfac));
            balls[2]._order=3;
            //brown.
            balls[3]._shape.setFillColor(sf::Color(131,87,43));
            balls[3]._x=brown_x;
            balls[3]._y=brown_y;
            spots[2].setPosition(sf::Vector2f(brown_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-brown_y*_sfac));
            balls[3]._order=4;
            //blue.
            balls[4]._shape.setFillColor(sf::Color(0,0,255));
            balls[4]._x=blue_x;
            balls[4]._y=blue_y;
            spots[3].setPosition(sf::Vector2f(blue_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-blue_y*_sfac));
            balls[4]._order=5;
            //pink.
            balls[5]._shape.setFillColor(sf::Color(255,105,180));
            balls[5]._x=pink_x;
            balls[5]._y=pink_y;
            spots[4].setPosition(sf::Vector2f(pink_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-pink_y*_sfac));
            balls[5]._order=6;
            //black.
            balls[6]._shape.setFillColor(sf::Color(0,0,0));
            balls[6]._x=black_x;
            balls[6]._y=black_y;
            spots[5].setPosition(sf::Vector2f(black_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-black_y*_sfac));
            balls[6]._order=7;

            double x;
            double y;
            int i=7;
            int gxmin;
            int gxmax;
            int gymin;
            int gymax;

            if (gametype!=3)
            {
                for (int row=0;row<5;row++)
                {
                    for (int h=-row;h<row+1;h=h+2)
                    {
                        balls[i]._shape.setFillColor(sf::Color(255,0,0));
                        balls[i]._x=pink_x-0.1-2*ball_radius-row*sqrt(3.0)*(ball_radius+DOUBLE_EPSILON);
                        balls[i]._y=pink_y+h*(ball_radius+DOUBLE_EPSILON);
                        balls[i]._order=i+1;
                        i+=1;
                    }
                }
            }
            else
            {
                for (int j=0;j<2;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=rail_thickness+cush_thickness+(black_x-rail_thickness-cush_thickness)*(j+1)/3.;
                    i+=1;
                }
                for (int j=0;j<5;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=black_x+(pink_x-black_x)*(j+1)/6.;
                    i+=1;
                }
                for (int j=0;j<5;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=pink_x+(blue_x-pink_x)*(j+1)/6.;
                    i+=1;
                }
                for (int j=0;j<3;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=blue_x+(blue_x-pink_x)*(j+1)/6.;
                    i+=1;
                }
            }

            for (int i=0;i<22;i++)
            {
                gxmin=int(floor((balls[i]._x-balls[i]._r)/(2.*balls[i]._r)));
                gxmax=int(floor((balls[i]._x+balls[i]._r)/(2.*balls[i]._r)));
                gymin=int(floor((balls[i]._y-balls[i]._r)/(2.*balls[i]._r)));
                gymax=int(floor((balls[i]._y+balls[i]._r)/(2.*balls[i]._r)));

                balls[i]._gpos[0][0]=gxmin;
                balls[i]._gpos[0][1]=gymin;
                balls[i]._gpos[1][0]=gxmax;
                balls[i]._gpos[1][1]=gymin;
                balls[i]._gpos[2][0]=gxmax;
                balls[i]._gpos[2][1]=gymax;
                balls[i]._gpos[3][0]=gxmin;
                balls[i]._gpos[3][1]=gymax;
            }

            baulkline.setPrimitiveType(sf::PrimitiveType::LineStrip);
            baulkline.resize(2);
            baulkcircle.setPrimitiveType(sf::PrimitiveType::LineStrip);
            baulkcircle.resize(100);
            baulkline[0].position=sf::Vector2f(_sfac*brown_x,_sfac*(rail_thickness+cush_thickness));
            baulkline[1].position=sf::Vector2f(_sfac*brown_x,(_sfac*raw_height/(1.+panel_ratio))-_sfac*(rail_thickness+cush_thickness));

            for (int i=0;i<100;i++)
            {
                x=brown_x+11.687*sin(pi*i/99);
                y=brown_y+11.687*cos(pi*i/99);
                baulkcircle[i].position=sf::Vector2f(_sfac*x,_sfac*y);
            }

            _shapes.push_back(&baulkline);
            _shapes.push_back(&baulkcircle);

            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&spots[i]);
            }

            cuetraj.setPrimitiveType(sf::PrimitiveType::LineStrip);
            cuetraj2.setPrimitiveType(sf::PrimitiveType::LineStrip);
            obtraj.setPrimitiveType(sf::PrimitiveType::LineStrip);
            ghostball.setRadius(ball_radius*_sfac);
            ghostball.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            ghostball.setFillColor(sf::Color(255,255,255,100));
            ghostball.setPosition(sf::Vector2f(-100.,-100.));
            ghostball2.setRadius(ball_radius*_sfac);
            ghostball2.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            ghostball2.setFillColor(sf::Color(255,255,255,100));
            ghostball2.setPosition(sf::Vector2f(-100.,-100.));

            _shapes.push_back(&cuetraj);
            _shapes.push_back(&cuetraj2);
            _shapes.push_back(&obtraj);
            _shapes.push_back(&ghostball);
            _shapes.push_back(&ghostball2);

            for (int i=0;i<22;i++)
            {
                _shapes.push_back(&balls[i]._shape);
            }

            _shapes.push_back(&cue._sprite);
            _shapes.push_back(&pauserect);
            _shapes.push_back(&pausetext);

            _shapes.push_back(&gameoverrect);
            _shapes.push_back(&gameovertext);
            _shapes.push_back(&stats_title);
            for (int i=0;i<stats_text.size();i++)
            {
                _shapes.push_back(&stats_text[i]);
            }
        }
        void update(double dt);
};

void GameScreen::update(double dt)
{
    if (!gameover)
    {
        //listen for packets.
        if (gametype<2)
        {
            packet.clear();
            if (socket.receive(packet)==sf::Socket::Done)
            {
                packet >> packetId;

                if (packetId==0)
                {
                    //display cue trajectory prediction.
                }
                else if (packetId==1)
                {
                    packet >> resultsize;
                    result.clear();
                    for (int i=0;i<resultsize;i++)
                    {
                        for (int j=0;j<66;j++)
                        {
                            packet >> temp[j];
                        }
                        result.push_back(temp);
                    }
                    t=0;
                    done=false;
                }
                else if (packetId==2)
                {
                    //whose turn it is.
                }
            }
        }

        if (int(floor(t*100./framerate))<result.size())
        {
            dt=t*100./framerate-int(floor(t*100./framerate));
            for (int i=0;i<22;i++)
            {
                if (int(floor(t*100./framerate))==int(result.size()-1))
                {
                    //dont interpolate.
                    balls[i]._x=result[result.size()-1][i*3];
                    balls[i]._y=result[result.size()-1][i*3+1];
                    balls[i]._z=result[result.size()-1][i*3+2];
                }
                else
                {
                    balls[i]._x=result[int(floor(t*100./framerate))][i*3]+dt*(result[int(floor(t*100./framerate))+1][i*3]-result[int(floor(t*100./framerate))][i*3]);
                    balls[i]._y=result[int(floor(t*100./framerate))][i*3+1]+dt*(result[int(floor(t*100./framerate))+1][i*3+1]-result[int(floor(t*100./framerate))][i*3+1]);
                    balls[i]._z=result[int(floor(t*100./framerate))][i*3+2]+dt*(result[int(floor(t*100./framerate))+1][i*3+2]-result[int(floor(t*100./framerate))][i*3+2]);
                }
            }
        }
        if (!done)
        {
            t=t+1;
            if (int(floor(t*100./framerate))>=result.size())
            {
                done=true;
                change=true;
                for (int i=0;i<22;i++)
                {
                    balls[i]._x=result[result.size()-1][i*3];
                    balls[i]._y=result[result.size()-1][i*3+1];
                    balls[i]._z=result[result.size()-1][i*3+2];
                }
            }
        }

        if (ispaused)
        {
            sf::Vector2f pos=pauserect.getPosition();
            double ds=sqrt(fabs(pos.y)*2.*pausea)*dt;
            pauserect.setPosition(sf::Vector2f(0.,fmin(0.,pos.y+ds)));
            pos=_buttons[0]._shape.getPosition();
            _buttons[0]._shape.setPosition(sf::Vector2f(pos.x,fmin(0.5*_sfac*raw_height,pos.y+ds)));
            _buttons[0]._text.setPosition(sf::Vector2f(pos.x,fmin(0.5*_sfac*raw_height,pos.y+ds)));
            pos=pausetext.getPosition();
            pausetext.setPosition(sf::Vector2f(pos.x,fmin(0.15*_sfac*raw_height,pos.y+ds)));
        }
        else
        {
            sf::Vector2f pos=pauserect.getPosition();
            double ds=-sqrt(fabs(pos.y+_sfac*raw_height)*2.*pausea)*dt;
            pauserect.setPosition(sf::Vector2f(0.,fmax(-_sfac*raw_height,pos.y+ds)));
            pos=_buttons[0]._shape.getPosition();
            _buttons[0]._shape.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            _buttons[0]._text.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            pos=pausetext.getPosition();
            pausetext.setPosition(sf::Vector2f(pos.x,fmax(-0.85*_sfac*raw_height,pos.y+ds)));
        }

        if (placing_white)
        {
            touching=false;
            for (int i=1;i<22;i++)
            {
                if (balls[0]._potted==false)
                {
                    if (sqrt(pow(balls[0]._x-balls[i]._x,2)+pow(balls[0]._y-balls[i]._y,2))<2*ball_radius)
                    {
                        balls[0]._shape.setFillColor(sf::Color(200,0,0));
                        touching=true;
                        break;
                    }
                }
            }
            if (touching==false)
            {
                balls[0]._shape.setFillColor(sf::Color(255,255,255));
            }
        }

        //get inputs.
        if (sf::Keyboard::isKeyPressed(user_controls["Pause game"]))
        {
            if (!ispausepressed)
            {
                ispausepressed=true;
                ispaused=!ispaused;
            }
        }
        else
        {
            ispausepressed=false;
        }

        change=false;
        if (isyourturn && done && fabs(pauserect.getPosition().y+_sfac*raw_height)<0.001)
        {
            //get user input when it is their turn.
            if (placing_white)
            {
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball up"]))
                {
                    //move ball up.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y+0.1*(dt/0.01)-brown_y,2.))<11.687)
                        {
                            balls[0]._y+=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._y=std::min(balls[0]._y+0.1*(dt/0.01),rail_thickness+cush_thickness+table_width-ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball down"]))
                {
                    //move ball down.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y-0.1*(dt/0.01)-brown_y,2.))<11.687)
                        {
                            balls[0]._y-=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._y=std::max(balls[0]._y-0.1*(dt/0.01),rail_thickness+cush_thickness+ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball left"]))
                {
                    //move ball left.
                    change=true;
                    if (gametype!=3)
                    {
                        balls[0]._x=std::max(balls[0]._x-0.1*(dt/0.01),brown_x);
                    }
                    else
                    {
                        balls[0]._x=std::max(balls[0]._x-0.1*(dt/0.01),rail_thickness+cush_thickness+ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball right"]))
                {
                    //move ball right.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x+0.1*(dt/0.01)-brown_x,2.)+pow(balls[0]._y-brown_y,2.))<11.687)
                        {
                            balls[0]._x+=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._x=std::min(balls[0]._x+0.1*(dt/0.01),rail_thickness+cush_thickness+table_length-ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Place ball"]))
                {
                    //place ball.
                    change=true;
                    if (!touching)
                    {
                        placing_white=false;
                    }
                }
            }
            else
            {
                if (sf::Keyboard::isKeyPressed(user_controls["Aim left"]))
                {
                    //move aim left.
                    change=true;
                    cue._angle-=(pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Aim right"]))
                {
                    //move aim right.
                    change=true;
                    cue._angle+=(pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Precise aim left"]))
                {
                    //precise aim left.
                    change=true;
                    cue._angle-=(0.015*pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Precise aim right"]))
                {
                    //precise aim right.
                    change=true;
                    cue._angle+=(0.015*pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Increase power"]))
                {
                    //increase power.
                    change=true;
                    power=fmin(power+0.5*(dt/0.01),100.5);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Decrease power"]))
                {
                    //decrease power.
                    change=true;
                    power=fmax(power-0.5*(dt/0.01),0.);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Increase cue elevation"]))
                {
                    //increase cue elevation.
                    change=true;
                    cue._alpha=fmin(cue._alpha+(0.1*pi/180.)*(dt/0.01),0.5*pi);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Decrease cue elevation"]))
                {
                    //decrease cue elevation.
                    change=true;
                    cue._alpha=fmax(cue._alpha-(0.1*pi/180.)*(dt/0.01),0.);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip up"]))
                {
                    //offset cue tip up.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta),2.)+pow(cue._offset*cos(cue._theta)+dt*ball_radius,2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta),cue._offset*cos(cue._theta)+dt*ball_radius);
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip down"]))
                {
                    //offset cue tip down.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta),2.)+pow(cue._offset*cos(cue._theta)-dt*ball_radius,2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta),cue._offset*cos(cue._theta)-dt*ball_radius);
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip left"]))
                {
                    //offset cue tip left.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta)-dt*ball_radius,2.)+pow(cue._offset*cos(cue._theta),2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta)-dt*ball_radius,cue._offset*cos(cue._theta));
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip right"]))
                {
                    //offset cue tip right.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta)+dt*ball_radius,2.)+pow(cue._offset*cos(cue._theta),2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta)+dt*ball_radius,cue._offset*cos(cue._theta));
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Strike cueball"]))
                {
                    //strike cueball.
                }
            }
        }

        for (int i=0;i<22;i++)
        {
            balls[i]._shape.setPosition(sf::Vector2f(_sfac*balls[i]._x,_sfac*raw_height/(1.+panel_ratio)-_sfac*balls[i]._y));
        }
        if (change && !placing_white)
        {
            sf::Color c;
            for (int i=0;i<100;i++)
            {
                c=power_bar[99-i].getFillColor();
                if (i<int(floor(power))) {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,255));}
                else {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,0));}
            }

            spin_selector_pos=spin_selector.getPosition();
            spin_dot.setPosition(sf::Vector2f(spin_selector_pos.x+cue._offset*sin(cue._theta)*spin_selector.getRadius()/ball_radius,spin_selector_pos.y-cue._offset*cos(cue._theta)*spin_selector.getRadius()/ball_radius));

            cue._sprite.setPosition(sf::Vector2f((balls[0]._x-2.*sin(cue._angle))*dfactor,window_height-dfactor*(balls[0]._y-2.*cos(cue._angle))));
            cue._sprite.setRotation((cue._angle+0.5*pi)*180./pi);

            cue._speed=1.2*power;
            cue.shot();
            balls[0]._vx=cue._ballv*sin(cue._angle);
            balls[0]._vy=cue._ballv*cos(cue._angle);
            balls[0]._xspin=cue._ballparspin*sin(cue._angle)+cue._ballperspin*cos(cue._angle);
            balls[0]._yspin=cue._ballparspin*cos(cue._angle)-cue._ballperspin*sin(cue._angle);
            balls[0]._rspin=cue._ballrspin;

            predict=trajectory(balls,cushions);

            cuetraj.clear();
            cuetraj2.clear();
            obtraj.clear();

            cuetraj.resize(predict[0].size()/2);
            cuetraj2.resize(1+predict[1].size()/2);
            obtraj.resize(predict[2].size()/2);

            if (predict[0].size()>1)
            {
                ghostball.setPosition(sf::Vector2f(dfactor*predict[0][predict[0].size()-2],window_height-dfactor*predict[0][predict[0].size()-1]));
                cuetraj2[0].position=sf::Vector2f(dfactor*predict[0][predict[0].size()-2],window_height-dfactor*predict[0][predict[0].size()-1]);
            }
            else
            {
                ghostball.setPosition(sf::Vector2f(-100.,-100.));
            }

            if (predict[1].size()>1)
            {
                ghostball2.setPosition(sf::Vector2f(dfactor*predict[1][predict[1].size()-2],window_height-dfactor*predict[1][predict[1].size()-1]));
            }
            else
            {
                ghostball2.setPosition(sf::Vector2f(-100.,-100.));
            }

            for (int i=0;i<predict[0].size()/2;i++)
            {
                cuetraj[i].position=sf::Vector2f(dfactor*predict[0][2*i],window_height-dfactor*predict[0][2*i+1]);
            }
            for (int i=1;i<1+predict[1].size()/2;i++)
            {
                cuetraj2[i].position=sf::Vector2f(dfactor*predict[1][2*i-2],window_height-dfactor*predict[1][2*i-1]);
            }
            for (int i=0;i<predict[2].size()/2;i++)
            {
                obtraj[i].position=sf::Vector2f(dfactor*predict[2][2*i],window_height-dfactor*predict[2][2*i+1]);
            }
        }
        if (placing_white || !done)
        {
            cue._sprite.setPosition(sf::Vector2f(-1000.,-1000.));
            ghostball.setPosition(sf::Vector2f(-1000.,-1000.));
            ghostball2.setPosition(sf::Vector2f(-1000.,-1000.));
            cuetraj.clear();
            cuetraj2.clear();
            obtraj.clear();
        }
    }
    else
    {
        //game over!
        sf::Vector2f pos=gameoverrect.getPosition();
        double ds=sqrt(fabs(pos.y)*2.*pausea)*dt;
        gameoverrect.setPosition(sf::Vector2f(0.,fmin(0.,pos.y+ds)));
        pos=_buttons[1]._shape.getPosition();
        _buttons[1]._shape.setPosition(sf::Vector2f(pos.x,fmin(0.75*_sfac*raw_height,pos.y+ds)));
        _buttons[1]._text.setPosition(sf::Vector2f(pos.x,fmin(0.75*_sfac*raw_height,pos.y+ds)));
        pos=gameovertext.getPosition();
        gameovertext.setPosition(sf::Vector2f(pos.x,fmin(0.15*_sfac*raw_height,pos.y+ds)));

        for (int i=0;i<stats_text.size();i++)
        {
            pos=stats_text[i].getPosition();
            stats_text[i].setPosition(sf::Vector2f(pos.x,fmin(stats_text_y[i]+_sfac*raw_height,pos.y+ds)));
        }
    }
}

#endif // CLIENT_H_INCLUDED
