#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <koolplot.h>
#include <string>
#include "objects.h"

int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    sf::ContextSettings settings;
    settings.antialiasingLevel=8;
    sf::RenderWindow window(sf::VideoMode(window_width,window_height+panel_height),"Snooker Game",sf::Style::Default,settings);

    //shaders.
    sf::Shader shader;
    if (!shader.loadFromFile("vertex_shader.vert","fragment_shader.frag"))
    {
        std::cout << "Error. Shader not loaded." << std::endl;
    }

    //set up panel.
    sf::RectangleShape panel_background(sf::Vector2f(window_width,panel_height));
    panel_background.setPosition(sf::Vector2f(0.,window_height));
    panel_background.setFillColor(sf::Color(0,0,0));

    sf::CircleShape spin_selector;
    spin_selector.setRadius(panel_height*0.4);
    spin_selector.setOrigin(panel_height*0.4,panel_height*0.4);
    spin_selector.setPosition(sf::Vector2f(0.5*panel_height,window_height+0.5*panel_height));
    spin_selector.setFillColor(sf::Color(255,255,255));
    spin_selector.setPointCount(200);

    sf::CircleShape spin_dot;
    spin_dot.setRadius(panel_height*0.4*0.005/0.02625);
    spin_dot.setOrigin(sf::Vector2f(panel_height*0.4*0.005/0.02625,panel_height*0.4*0.005/0.02625));
    spin_dot.setPosition(sf::Vector2f(0.5*panel_height,window_height+0.5*panel_height));
    spin_dot.setFillColor(sf::Color(255,0,0));
    spin_dot.setPointCount(100);

    std::array<sf::RectangleShape,100> power_bar;
    for (int i=0;i<100;i++)
    {
        power_bar[i].setSize(sf::Vector2f(0.1*panel_height,0.008*panel_height));
        power_bar[i].setPosition(sf::Vector2f(1.1*panel_height,window_height+0.1*panel_height+i*0.008*panel_height));
        power_bar[i].setFillColor(sf::Color(255,int(floor(2.55*i)),0));
    }
    double outline=0.01*panel_height;
    sf::RectangleShape power_outline(sf::Vector2f(0.1*panel_height+2*outline,0.8*panel_height+2*outline));
    power_outline.setPosition(sf::Vector2f(1.1*panel_height-outline,window_height+0.1*panel_height-outline));
    power_outline.setFillColor(sf::Color(0,0,0));
    power_outline.setOutlineThickness(outline);
    power_outline.setOutlineColor(sf::Color(255,255,255));

    sf::CircleShape elevation_border;
    elevation_border.setRadius(panel_height*0.4);
    elevation_border.setOrigin(panel_height*0.4,panel_height*0.4);
    elevation_border.setPosition(sf::Vector2f(window_width-0.5*panel_height,window_height+0.5*panel_height));
    elevation_border.setFillColor(sf::Color(0,0,0));
    elevation_border.setOutlineThickness(outline);
    elevation_border.setOutlineColor(sf::Color(255,255,255));
    elevation_border.setPointCount(200);

    sf::CircleShape elevation_ball;
    elevation_ball.setRadius(0.2*panel_height);
    elevation_ball.setOrigin(0.2*panel_height,0.2*panel_height);
    elevation_ball.setPosition(sf::Vector2f(window_width-0.5*panel_height,window_height+0.5*panel_height));
    elevation_ball.setFillColor(sf::Color(255,255,255));
    elevation_ball.setPointCount(150);

    sf::RectangleShape elevation_pointer(sf::Vector2f(0.2*panel_height,0.2*panel_height*0.005/0.02625));
    elevation_pointer.setOrigin(0,0.1*panel_height*0.005/0.02625);
    elevation_pointer.setPosition(sf::Vector2f(window_width-0.3*panel_height,window_height+0.5*panel_height));
    elevation_pointer.setFillColor(sf::Color(255,0,0));

    sf::Font font;
    if (!font.loadFromFile("Roboto-Thin.ttf"))
    {
        std::cout << "ERROR. Couldn't load font!" << std::endl;
    }
    sf::Text elevation_display;
    elevation_display.setFont(font);
    elevation_display.setString("00�");
    elevation_display.setCharacterSize(int(0.15*panel_height));
    elevation_display.setFillColor(sf::Color(0,0,0));
    sf::Vector2f pos=elevation_display.getScale();
    elevation_display.setOrigin(sf::Vector2f(0.5*0.15*panel_height*pos.x/pos.y,0.5*0.15*panel_height));
    elevation_display.setPosition(sf::Vector2f(window_width-0.5*panel_height,window_height+0.5*panel_height));

    //get cue.
    Cue cue;
    cue._sprite.setPosition(sf::Vector2f(black_x*dfactor,window_height-black_y*dfactor));

    Ball balls[22];

    //prepare cushions.
    Cushion cushions[6];

    cushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
    cushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
    cushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
    cushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
    cushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
    cushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

    //prepare matrix.
    for (int i=0;i<44;i++)
    {
        M_(i,i)=1/ball_mass;
    }
    M_(44,44)=0;
    M_(45,45)=0;

    //prepare the balls and spots.

    sf::CircleShape spots[6];

    for (int i=0;i<6;i++)
    {
        spots[i].setRadius(spot_r*dfactor);
        spots[i].setOrigin(spot_r*dfactor,spot_r*dfactor);
        spots[i].setFillColor(sf::Color(255,255,255));
    }

    for (int i=0;i<22;i++)
    {
        balls[i]._shape.setOrigin(ball_radius*dfactor,ball_radius*dfactor);
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
    spots[0].setPosition(sf::Vector2f(yellow_x*dfactor,window_height-yellow_y*dfactor));
    balls[1]._order=2;
    //green.
    balls[2]._shape.setFillColor(sf::Color(0,150,0));
    balls[2]._x=green_x;
    balls[2]._y=green_y;
    spots[1].setPosition(sf::Vector2f(green_x*dfactor,window_height-green_y*dfactor));
    balls[2]._order=3;
    //brown.
    balls[3]._shape.setFillColor(sf::Color(131,87,43));
    balls[3]._x=brown_x;
    balls[3]._y=brown_y;
    spots[2].setPosition(sf::Vector2f(brown_x*dfactor,window_height-brown_y*dfactor));
    balls[3]._order=4;
    //blue.
    balls[4]._shape.setFillColor(sf::Color(0,0,255));
    balls[4]._x=blue_x;
    balls[4]._y=blue_y;
    spots[3].setPosition(sf::Vector2f(blue_x*dfactor,window_height-blue_y*dfactor));
    balls[4]._order=5;
    //pink.
    balls[5]._shape.setFillColor(sf::Color(255,105,180));
    balls[5]._x=pink_x;
    balls[5]._y=pink_y;
    spots[4].setPosition(sf::Vector2f(pink_x*dfactor,window_height-pink_y*dfactor));
    balls[5]._order=6;
    //black.
    balls[6]._shape.setFillColor(sf::Color(0,0,0));
    balls[6]._x=black_x;
    balls[6]._y=black_y;
    spots[5].setPosition(sf::Vector2f(black_x*dfactor,window_height-black_y*dfactor));
    balls[6]._order=7;

    double x;
    double y;
    int i=7;
    int gxmin;
    int gxmax;
    int gymin;
    int gymax;
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

    //draw lines.
    sf::Vertex baulkline[2];
    sf::VertexArray baulkcircle(sf::PrimitiveType::LineStrip,100);

    baulkline[0]=sf::Vertex(sf::Vector2f(dfactor*brown_x,dfactor*(rail_thickness+cush_thickness)));
    baulkline[1]=sf::Vertex(sf::Vector2f(dfactor*brown_x,window_height-dfactor*(rail_thickness+cush_thickness)));

    for (int i=0;i<100;i++)
    {
        x=brown_x+11.687*sin(pi*i/99);
        y=brown_y+11.687*cos(pi*i/99);
        baulkcircle[i].position=sf::Vector2f(dfactor*x,dfactor*y);
    }

    //final variables.
    Eigen::MatrixXd test=Eigen::MatrixXd::Constant(46,1,0.0);

    bool placing_white=false;
    bool touching=false;
    double power=50.;

    sf::Vector2i mouse_pos;
    sf::Vector2f spin_selector_pos;
    sf::Vector2f spin_dot_pos;
    double dist;
    double angle;
    double alpha=0.;
    std::array<double,2> roots;

    sf::Time elapsed;
    sf::Time diff;
    sf::Time period=sf::seconds(1.0/framerate);
    sf::Clock clock;

    std::vector<std::array<double,66>> result;
    result=simulate(balls,cushions);

    int t=0;

    while (window.isOpen())
    {
        if (t<result.size())
        {
            for (int i=0;i<22;i++)
            {
                balls[i]._x=result[t][i*3];
                balls[i]._y=result[t][i*3+1];
                balls[i]._z=result[t][i*3+2];
            }
            if (t==int(result.size()-1))
            {
                for (int i=0;i<22;i++)
                {
                    balls[i]._rspin=0;
                    if (balls[i]._x<0 && i<7 && i>0)
                    {
                        //respot colour.
                        balls[i]._potted=false;
                        balls[i]._z=ball_radius;
                        balls[i]._vx=0.;
                        balls[i]._vy=0.;
                        balls[i]._vz=0.;
                        balls[i]._xspin=0.;
                        balls[i]._yspin=0.;

                        if (i==1)
                        {
                            balls[i]._x=yellow_x;
                            balls[i]._y=yellow_y;
                        }
                        else if (i==2)
                        {
                            balls[i]._x=green_x;
                            balls[i]._y=green_y;
                        }
                        else if (i==3)
                        {
                            balls[i]._x=brown_x;
                            balls[i]._y=brown_y;
                        }
                        else if (i==4)
                        {
                            balls[i]._x=blue_x;
                            balls[i]._y=blue_y;
                        }
                        else if (i==5)
                        {
                            balls[i]._x=pink_x;
                            balls[i]._y=pink_y;
                        }
                        else if (i==6)
                        {
                            balls[i]._x=black_x;
                            balls[i]._y=black_y;
                        }
                    }
                }
            }
        }
        t=std::min(t+1,int(result.size()));

        //check if white ball potted.
        if (t==int(result.size()) && balls[0]._x<0)
        {
            placing_white=true;
            balls[0]._potted=false;
            balls[0]._x=cueball_break_x;
            balls[0]._y=cueball_break_y;
            balls[0]._z=ball_radius;
            balls[0]._vx=0;
            balls[0]._vy=0;
            balls[0]._vz=0;
            balls[0]._rspin=0;
            balls[0]._xspin=0;
            balls[0]._yspin=0;
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

        //deal with user input.
        sf::Event event;
        while (window.pollEvent(event))
        {
            switch (event.type)
            {
                case sf::Event::Closed:
                    window.close();
                    break;
                case sf::Event::KeyPressed:
                    if (event.key.code==sf::Keyboard::Escape)
                    {
                        window.close();
                        break;
                    }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
        {
            //left arrow.
            if (placing_white)
            {
                balls[0]._x=std::max(balls[0]._x-0.1,brown_x);
            }
            else if (t==int(result.size()))
            {
                //moving cue.
                cue._angle-=pi/180.;
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
        {
            //right arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x+0.1-brown_x,2.)+pow(balls[0]._y-brown_y,2.))<11.687)
                {
                    balls[0]._x+=0.1;
                }
            }
            else if (t==int(result.size()))
            {
                cue._angle+=pi/180.;
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
        {
            //up arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y+0.1-brown_y,2.))<11.687)
                {
                    balls[0]._y+=0.1;
                }
            }
            else if (t==int(result.size()))
            {
                power=fmin(power+0.5,100.5);
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
        {
            //down arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y-0.1-brown_y,2.))<11.687)
                {
                    balls[0]._y-=0.1;
                }
            }
            else if (t==int(result.size()))
            {
                power=fmax(power-0.5,0.);
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Enter))
        {
            if (placing_white)
            {
                if (balls[0]._shape.getFillColor()==sf::Color(255,255,255))
                {
                    placing_white=false;
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space))
        {
            if (!placing_white && t==int(result.size()))
            {
                //take the shot!
                spin_selector_pos=spin_selector.getPosition();
                spin_dot_pos=spin_dot.getPosition();
                dist=(ball_radius/(0.4*panel_height))*sqrt(pow(spin_selector_pos.x-spin_dot_pos.x,2.)+pow(spin_selector_pos.y-spin_dot_pos.y,2.));
                angle=atan2(spin_dot_pos.x-spin_selector_pos.x,spin_selector_pos.y-spin_dot_pos.y);
                cue._offset=dist;
                cue._theta=angle;
                cue._speed=1.2*power;
                cue.shot();
                balls[0]._vx=cue._ballv*sin(cue._angle);
                balls[0]._vy=cue._ballv*cos(cue._angle);
                balls[0]._xspin=cue._ballparspin*sin(cue._angle)+cue._ballperspin*cos(cue._angle);
                balls[0]._yspin=cue._ballparspin*cos(cue._angle)-cue._ballperspin*sin(cue._angle);
                balls[0]._rspin=cue._ballrspin;

                result=simulate(balls,cushions);
                t=0;
                continue;
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Backspace))
        {
            //replay for debug purposes.
            t=0;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Comma))
        {
            if (!placing_white && t==int(result.size()))
            {
                //move left.
                cue._angle-=0.1*pi/180.;
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Period))
        {
            if (!placing_white && t==int(result.size()))
            {
                //move right.
                cue._angle+=0.1*pi/180.;
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::LBracket))
        {
            if (!placing_white && t==int(result.size()))
            {
                //elevate the cue upwards.
                alpha=fmin(alpha+0.1*pi/180.,0.5*pi);

                if (int(alpha*180./pi)<10)
                {
                    elevation_display.setString("0"+std::to_string(int(alpha*180./pi))+"�");
                }
                else
                {
                    elevation_display.setString(std::to_string(int(alpha*180./pi))+"�");
                }

                spin_selector_pos=spin_selector.getPosition();
                spin_dot_pos=spin_dot.getPosition();
                dist=(spin_selector_pos.y-spin_dot_pos.y)*0.5;
                roots=qsolve_quadratic(1.,2*(sqrt(pow(0.2*panel_height,2.)-pow(dist,2.))*cos(alpha)+dist*sin(alpha)),pow(0.2*panel_height,2.)-pow(0.4*panel_height,2.));
                elevation_pointer.setPosition(sf::Vector2f(window_width-0.5*panel_height+sqrt(pow(0.2*panel_height,2.)-pow(dist,2.)),window_height+0.5*panel_height-dist));
                elevation_pointer.setRotation(-alpha*180./pi);
                if (roots[0]!=roots[0])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[1],0.2*panel_height*0.005/0.02625));
                }
                else if (roots[1]!=roots[1])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[0],0.2*panel_height*0.005/0.02625));
                }
                else
                {
                    elevation_pointer.setSize(sf::Vector2f(fmax(roots[0],roots[1]),0.2*panel_height*0.005/0.02625));
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::RBracket))
        {
            if (!placing_white && t==int(result.size()))
            {
                //flatten the cue.
                alpha=fmax(alpha-0.1*pi/180.,0.);

                if (int(alpha*180./pi)<10)
                {
                    elevation_display.setString("0"+std::to_string(int(alpha*180./pi))+"�");
                }
                else
                {
                    elevation_display.setString(std::to_string(int(alpha*180./pi))+"�");
                }

                spin_selector_pos=spin_selector.getPosition();
                spin_dot_pos=spin_dot.getPosition();
                dist=(spin_selector_pos.y-spin_dot_pos.y)*0.5;
                roots=qsolve_quadratic(1.,2*(sqrt(pow(0.2*panel_height,2.)-pow(dist,2.))*cos(alpha)+dist*sin(alpha)),pow(0.2*panel_height,2.)-pow(0.4*panel_height,2.));
                elevation_pointer.setPosition(sf::Vector2f(window_width-0.5*panel_height+sqrt(pow(0.2*panel_height,2.)-pow(dist,2.)),window_height+0.5*panel_height-dist));
                elevation_pointer.setRotation(-alpha*180./pi);
                if (roots[0]!=roots[0])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[1],0.2*panel_height*0.005/0.02625));
                }
                else if (roots[1]!=roots[1])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[0],0.2*panel_height*0.005/0.02625));
                }
                else
                {
                    elevation_pointer.setSize(sf::Vector2f(fmax(roots[0],roots[1]),0.2*panel_height*0.005/0.02625));
                }
            }
        }
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
        {
            mouse_pos=sf::Mouse::getPosition(window);
            if (mouse_pos.x>=0 && mouse_pos.x<=panel_height && mouse_pos.y>=window_height && mouse_pos.y<=window_height+panel_height)
            {
                dist=fmin(sqrt(pow(mouse_pos.x-0.5*panel_height,2.)+pow(mouse_pos.y-0.5*panel_height-window_height,2.)),0.4*panel_height-panel_height*0.4*0.005/0.02625);
                angle=atan2(mouse_pos.x-0.5*panel_height,mouse_pos.y-0.5*panel_height-window_height);
                spin_dot.setPosition(sf::Vector2f(0.5*panel_height+dist*sin(angle),window_height+panel_height*0.5+dist*cos(angle)));

                spin_selector_pos=spin_selector.getPosition();
                spin_dot_pos=spin_dot.getPosition();
                dist=(spin_selector_pos.y-spin_dot_pos.y)*0.5;
                roots=qsolve_quadratic(1.,2*(sqrt(pow(0.2*panel_height,2.)-pow(dist,2.))*cos(alpha)+dist*sin(alpha)),pow(0.2*panel_height,2.)-pow(0.4*panel_height,2.));
                elevation_pointer.setPosition(sf::Vector2f(window_width-0.5*panel_height+sqrt(pow(0.2*panel_height,2.)-pow(dist,2.)),window_height+0.5*panel_height-dist));
                elevation_pointer.setRotation(-alpha*180./pi);
                if (roots[0]!=roots[0])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[1],0.2*panel_height*0.005/0.02625));
                }
                else if (roots[1]!=roots[1])
                {
                    elevation_pointer.setSize(sf::Vector2f(roots[0],0.2*panel_height*0.005/0.02625));
                }
                else
                {
                    elevation_pointer.setSize(sf::Vector2f(fmax(roots[0],roots[1]),0.2*panel_height*0.005/0.02625));
                }
            }
        }

        //draw all the shapes on the screen.
        window.clear(sf::Color(baizecolour[0],baizecolour[1],baizecolour[2]));

        window.draw(panel_background);
        window.draw(spin_selector);
        window.draw(power_outline);
        window.draw(spin_dot);
        window.draw(elevation_border);
        window.draw(elevation_ball);
        window.draw(elevation_pointer);
        window.draw(elevation_display);

        for (int i=0;i<int(floor(power));i++)
        {
            window.draw(power_bar[99-i]);
        }

        for (int i=0;i<6;i++)
        {
            window.draw(spots[i]);
            window.draw(cushions[i]._pocketshape);
        }
        window.draw(baulkline,2,sf::Lines);
        window.draw(baulkcircle);
        for (int i=0;i<6;i++)
        {
            window.draw(cushions[i]._shape);
            window.draw(cushions[i]._railshape);
            window.draw(cushions[i]._p1shape);
            window.draw(cushions[i]._p2shape);
        }
        for (int i=0;i<22;i++)
        {
            balls[i]._shape.setPosition(sf::Vector2f(balls[i]._x*dfactor,window_height-balls[i]._y*dfactor));
            window.draw(balls[i]._shape,&shader);
        }
        if (!placing_white && t==int(result.size()))
        {
            cue._sprite.setPosition(sf::Vector2f((balls[0]._x-2.*sin(cue._angle))*dfactor,window_height-dfactor*(balls[0]._y-2.*cos(cue._angle))));
            cue._sprite.setRotation((cue._angle+0.5*pi)*180./pi);
            window.draw(cue._sprite);
        }

        elapsed=clock.restart();
        diff=period-elapsed;
        if (diff.asSeconds()>0)
        {
            sf::sleep(diff);
        }

        window.display();
    }

    return 0;
}
