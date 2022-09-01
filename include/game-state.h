#ifndef GAME-STATE_H_INCLUDED
#define GAME-STATE_H_INCLUDED

class GameState
{
    public:
        bool _shouldUpdate=true;

        std::vector<RectButton> _buttons;
        std::vector<InputBox> _inputboxes;
        std::vector<sf::Drawable*> _shapes;
        std::vector<sf::Drawable*> _importantShapes;

        sf::Color _background;

        double _sfac;

        GameState(double sfac) {_sfac=sfac;}
        virtual void update(double dt,sf::Vector2i mouse_pos)=0;
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
