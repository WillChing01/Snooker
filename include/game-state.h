#ifndef GAME-STATE_H_INCLUDED
#define GAME-STATE_H_INCLUDED

class GameState
{
    public:
        bool _shouldUpdate=true;

        std::vector<RectButton> _buttons;
        std::vector<InputBox> _inputboxes;
        std::vector<sf::Drawable*> _shapes;

        sf::Color _background;

        double _sfac;

        GameState(double sfac) {_sfac=sfac;}
        virtual void update(double dt,sf::Vector2i mouse_pos)=0;
};

#endif // GAME-STATE_H_INCLUDED
