#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <koolplot.h>
#include "objects.h"

using namespace boost::numeric::odeint;

//int main()
//{
//    double v0=2.;
//    double delta=1.5;
//    double theta=asin(delta/mpocket_r);
//
//    std::vector<state_type> output;
//    std::vector<double> times;
//    state_type xi={0,v0*cos(theta)/ball_radius,theta,v0*sin(theta)/mpocket_r};
//
//    Rimode rimode(mpocket_r);
//
//    runge_kutta4<state_type> stepper;
//
//    size_t steps=integrate(rimode,xi,0.,1.,0.00001,push_back_state_and_time(output,times));
//    //integrate_const(stepper,rimode,xi,0.,1.,pow(10,-6),push_back_state_and_time(output,times));
//
//    std::cout << "Steps: " << steps << std::endl;
//
//    double normal;
//
//
//    for (int i=0;i<steps;i++)
//    {
//        normal=gravity*cos(output[i][0])-ball_radius*pow(output[i][1],2)+(mpocket_r-ball_radius*sin(output[i][0]))*sin(output[i][0])*pow(output[i][3],2);
//        if (normal<0)
//        {
//            break;
//        }
//        std::cout << times[i] << " : " << (mpocket_r-ball_radius*sin(output[i][0]))*sin(output[i][2]) << " : " << -(mpocket_r-ball_radius*sin(output[i][0]))*cos(output[i][2]) << std::endl;
//    }
//
//    return 0;
//}

int main()
{
    sf::ContextSettings settings;
    settings.antialiasingLevel=8;
    sf::RenderWindow window(sf::VideoMode(window_width,window_height),"Snooker Game",sf::Style::Default,settings);

    //shaders.
    sf::Shader shader;
    if (!shader.loadFromFile("vertex_shader.vert","fragment_shader.frag"))
    {
        std::cout << "Error. Shader not loaded." << std::endl;
    }

    //prepare the balls and spots.
    Ball balls[22];
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
    balls[0]._y=cueball_break_y-0.7;
    //balls[0]._y=rail_thickness+cush_thickness+ball_radius;
    //__float128 angle=atan2q(cpockets[0][0]-balls[0]._x,cpockets[0][1]-balls[0]._y)+0.5*pi*0.009;
    double angle=-0.5*pi+0.001;
    balls[0]._vx=75*sin(angle);
    balls[0]._vy=75*cos(angle);
    balls[0]._xspin=0;
    balls[0]._yspin=0;
    balls[0]._rspin=-80;
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
            balls[i]._x=pink_x-0.1-2*ball_radius-row*sqrt(3.0);
            balls[i]._y=pink_y+h*ball_radius;
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

    //prepare cushions.
    Cushion cushions[6];

    cushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
    cushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
    cushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
    cushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
    cushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
    cushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

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

    //prepare matrix.
    for (int i=0;i<44;i++)
    {
        M_(i,i)=1/ball_mass;
    }
    M_(44,44)=0;
    M_(45,45)=0;

    //final variables.
    double newx;
    double newy;
    Eigen::MatrixXd test=Eigen::MatrixXd::Constant(46,1,0.0);

    bool still=true;
    bool placing_white=false;
    bool touching=false;

    sf::Time elapsed;
    sf::Time diff;
    sf::Time period=sf::seconds(1.0/framerate);
    sf::Clock clock;

    std::vector<std::array<double,66>> result;
    result=simulate(balls,cushions);

    int t=0;

    while (window.isOpen())
    {
        for (int i=0;i<22;i++)
        {
            if (t<result.size())
            {
                balls[i]._x=result[t][i*3];
                balls[i]._y=result[t][i*3+1];
                balls[i]._z=result[t][i*3+2];
            }
            else
            {
                balls[i]._x=result[result.size()-1][i*3];
                balls[i]._y=result[result.size()-1][i*3+1];
                balls[i]._z=result[result.size()-1][i*3+2];
            }
        }
        t+=1;
        //check if white ball potted.
        if (balls[0]._potted && still)
        {
            placing_white=true;
            balls[0]._potted=false;
            balls[0]._x=cueball_break_x;
            balls[0]._y=cueball_break_y;
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
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
        {
            //right arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x+0.1-brown_x,2)+pow(balls[0]._y-brown_y,2))<11.687)
                {
                    balls[0]._x+=0.1;
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
        {
            //up arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x-brown_x,2)+pow(balls[0]._y+0.1-brown_y,2))<11.687)
                {
                    balls[0]._y+=0.1;
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
        {
            //down arrow.
            if (placing_white)
            {
                if (sqrt(pow(balls[0]._x-brown_x,2)+pow(balls[0]._y-0.1-brown_y,2))<11.687)
                {
                    balls[0]._y-=0.1;
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space))
        {
            if (placing_white)
            {
                if (balls[0]._shape.getFillColor()==sf::Color(255,255,255))
                {
                    placing_white=false;
                }
            }
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Backspace))
        {
            //replay for debug purposes.
            t=0;
        }

        //draw all the shapes on the screen.
        window.clear(sf::Color(baizecolour[0],baizecolour[1],baizecolour[2]));
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
            if (balls[i]._potted==false)
            {
                balls[i]._shape.setPosition(sf::Vector2f(balls[i]._x*dfactor,window_height-balls[i]._y*dfactor));
                window.draw(balls[i]._shape,&shader);
            }
        }

        elapsed=clock.restart();
        diff=period-elapsed;
        if (diff.asSeconds()>0)
        {
            sf::sleep(diff);
        }
        else
        {
            std::cout << "Lagging! Behind by:" << -diff.asSeconds() << std::endl;
        }

        window.display();
    }

    return 0;
}
