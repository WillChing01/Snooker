#ifndef GAME-STATE_H_INCLUDED
#define GAME-STATE_H_INCLUDED

const std::vector<std::string> _stateTargets=
{
    "Quit",
    "TitleScreen",
    "Singleplayer",
    "SingleplayerAI",
    "SingleplayerLineup",
    "Multiplayer",
    "MultiplayerHost",
    "MultiplayerJoin",
    "Options",
    "Controls",
    "Changecue"
};

class GameState
{
    public:
        bool _shouldUpdate=true;
        bool _shouldChangeState=false;
        std::string _stateTarget="";
        bool _isWaitingForInput=false;
        int _controlindex=-1;
        bool _isTyping=false;
        int _typingIndex=0;

        std::vector<RectButton> _buttons;
        std::vector<InputBox> _inputboxes;
        std::vector<sf::Drawable*> _shapes;
        std::vector<sf::Drawable*> _importantShapes;

        sf::Color _background;

        double _sfac;

        GameState(double sfac) {_sfac=sfac;}
        virtual void update(double dt,sf::Vector2i mouse_pos)=0;

        virtual void handleButtonCallbacks()
        {
            for (int i=0;i<_buttons.size();i++)
            {
                if (!_buttons[i]._isactive) {continue;}
                if (!_buttons[i]._shouldExecute) {continue;}

                _buttons[i]._shouldExecute=false;

                if (_buttons[i]._callback)
                {
                    _buttons[i]._callback(static_cast<GameState*>(this),&_buttons[i]._payload);
                }

                break;
            }
        }

        virtual void handleButtonPress(sf::RenderWindow &window)
        {
            //check if any buttons are pressed.
            double sdiff=sf::VideoMode::getDesktopMode().width-dfactor*raw_width;

            sf::Vector2i mouse_pos=sf::Mouse::getPosition(window);
            mouse_pos.x=int(mouse_pos.x-0.5*sdiff);
            bool isPressed=sf::Mouse::isButtonPressed(sf::Mouse::Left);

            if (!(mouse_pos.x>=0 && mouse_pos.x<=dfactor*raw_width && mouse_pos.y>=0 && mouse_pos.y<=dfactor*raw_height))
            {
                //mouse is outside of bounds.
                return;
            }

            sf::FloatRect bounds;
            for (int i=0;i<_buttons.size();i++)
            {
                //update the pressing state of the button.
                _buttons[i]._wasClicked=_buttons[i]._isClicked;
                _buttons[i]._isClicked=false;

                bounds=_buttons[i]._shape.getGlobalBounds();

                if (!bounds.contains(sf::Vector2f(mouse_pos)) || !_buttons[i]._isactive)
                {
                    _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                    _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                    _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
                    continue;
                }

                //now mouse is inside an active button.
                _buttons[i]._isClicked=isPressed;

                //set the appropriate colours.
                if (!isPressed)
                {
                    _buttons[i]._shape.setFillColor(_buttons[i]._colour2);
                    _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour2);
                    _buttons[i]._text.setFillColor(_buttons[i]._textcolour2);
                }
                else
                {
                    _buttons[i]._shape.setFillColor(_buttons[i]._colour3);
                    _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour3);
                    _buttons[i]._text.setFillColor(_buttons[i]._textcolour3);
                }

                //should update the state if a button was pressed.
                if (_buttons[i]._wasClicked && !_buttons[i]._isClicked)
                {
                    _buttons[i]._shouldExecute=true;
                    _shouldUpdate=true;
                    for (int j=0;j<_stateTargets.size();j++)
                    {
                        if (_buttons[i]._target==_stateTargets[j])
                        {
                            _shouldChangeState=true;
                            _stateTarget=_stateTargets[j];
                            _buttons[i]._shouldExecute=false; //not necessary for state change.
                            break;
                        }
                    }
                }
            }
        }

        virtual void handleInputClick(sf::RenderWindow &window)
        {
            double sdiff=sf::VideoMode::getDesktopMode().width-dfactor*raw_width;

            sf::Vector2i mouse_pos=sf::Mouse::getPosition(window);
            mouse_pos.x=int(mouse_pos.x-0.5*sdiff);
            bool isPressed=sf::Mouse::isButtonPressed(sf::Mouse::Left);

            if (!(mouse_pos.x>=0 && mouse_pos.x<=dfactor*raw_width && mouse_pos.y>=0 && mouse_pos.y<=dfactor*raw_height))
            {
                //mouse is outside of bounds.
                return;
            }

            sf::FloatRect bounds;
            for (int i=0;i<_inputboxes.size();i++)
            {
                bounds=_inputboxes[i]._shape.getGlobalBounds();

                if (!_inputboxes[i]._isactive || (isPressed && !bounds.contains(sf::Vector2f(mouse_pos))))
                {
                    //reset to default state.
                    _inputboxes[i]._isTyping=false;
                    _inputboxes[i]._shape.setOutlineColor(_inputboxes[i]._outlinecolour1);
                    sf::Color c=_inputboxes[i]._backtext.getFillColor();
                    if (_inputboxes[i]._input.empty())
                    {
                        _inputboxes[i]._backtext.setFillColor(sf::Color(c.r,c.g,c.b,150));
                    }
                    c=_inputboxes[i]._cursor.getFillColor();
                    _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,0));
                    continue;
                }

                if (isPressed && bounds.contains(sf::Vector2f(mouse_pos)))
                {
                    //set to typing state.
                    _inputboxes[i]._isTyping=true;
                    _inputboxes[i]._shape.setOutlineColor(_inputboxes[i]._outlinecolour2);
                    sf::Color c=_inputboxes[i]._backtext.getFillColor();
                    _inputboxes[i]._backtext.setFillColor(sf::Color(c.r,c.g,c.b,0));
                    c=_inputboxes[i]._cursor.getFillColor();
                    _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));

                    double dist=999.;
                    for (int j=0;j<_inputboxes[i]._input.length()+1;j++)
                    {
                        _inputboxes[i]._text.setString(_inputboxes[i]._input.substr(0,j));
                        bounds=_inputboxes[i]._text.getGlobalBounds();
                        if (fabs(mouse_pos.x-bounds.left-bounds.width)<dist)
                        {
                            dist=fabs(mouse_pos.x-bounds.left-bounds.width);
                            _inputboxes[i]._cursorpos=j;
                            sf::Vector2f rectpos=_inputboxes[i]._cursor.getPosition();
                            _inputboxes[i]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*_inputboxes[i]._abscursorthickness,rectpos.y));
                        }
                        _inputboxes[i]._text.setString(_inputboxes[i]._input);
                    }
                }
            }

            _isTyping=false;
            for (int i=0;i<_inputboxes.size();i++)
            {
                if (_inputboxes[i]._isTyping==true)
                {
                    _isTyping=true;
                    _typingIndex=i;
                    break;
                }
            }
        }

        virtual void render(sf::RenderWindow &window)
        {
            window.clear(_background);

            for (int i=0;i<_shapes.size();i++)
            {
                window.draw(*_shapes[i]);
            }
            for (int i=0;i<_buttons.size();i++)
            {
                window.draw(_buttons[i]._shape);
                window.draw(_buttons[i]._text);
            }
            for (int i=0;i<_inputboxes.size();i++)
            {
                window.draw(_inputboxes[i]._shape);
                window.draw(_inputboxes[i]._backtext);
                window.draw(_inputboxes[i]._text);
                window.draw(_inputboxes[i]._cursor);
            }
            for (int i=0;i<_importantShapes.size();i++)
            {
                window.draw(*_importantShapes[i]);
            }
            window.display();
        }
};

#endif // GAME-STATE_H_INCLUDED
