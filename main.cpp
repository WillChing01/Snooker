#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Network.hpp>
#include <string>
#include <thread>
#include "objects.h"
#include "client.h"

int main()
{
//    Server server;
//    std::thread server_thread(&Server::executionThread,std::ref(server));
//    server_thread.detach();
//
//    sf::TcpSocket socket;
//    if (socket.connect(sf::IpAddress(targetip),targetport)!=sf::Socket::Done)
//    {
//        //error.
//        std::cout << "Could not connect to self server" << std::endl;
//    }
//    socket.setBlocking(false);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    sf::ContextSettings settings;
    settings.antialiasingLevel=8;
    sf::RenderWindow window(sf::VideoMode(sf::VideoMode::getDesktopMode().width,sf::VideoMode::getDesktopMode().height),"Snooker Game",sf::Style::Fullscreen,settings);

    //shaders.
    sf::Shader shader;
    if (!shader.loadFromFile("vertex_shader.vert","fragment_shader.frag"))
    {
        std::cout << "Error. Shader not loaded." << std::endl;
    }

    double sdiff=sf::VideoMode::getDesktopMode().width-dfactor*raw_width;

    sf::RectangleShape p1(sf::Vector2f(0.5*sdiff,raw_height*dfactor));
    p1.setFillColor(sf::Color(0,0,0));
    p1.setPosition(sf::Vector2f(-0.5*sdiff,0.));

    sf::RectangleShape p2(sf::Vector2f(0.5*sdiff,raw_height*dfactor));
    p2.setFillColor(sf::Color(0,0,0));
    p2.setPosition(sf::Vector2f(dfactor*raw_width,0.));

    sf::View view;
    view.reset(sf::FloatRect(-0.5*sdiff,0.,float(sf::VideoMode::getDesktopMode().width),float(sf::VideoMode::getDesktopMode().height)));
    //view.setCenter(float(0.5*sf::VideoMode::getDesktopMode().width-0.5*sdiff),float(0.5*sf::VideoMode::getDesktopMode().height));
    window.setView(view);

    //prepare matrix.
    for (int i=0;i<44;i++)
    {
        M_(i,i)=1/ball_mass;
    }
    M_(44,44)=0;
    M_(45,45)=0;

    sf::Time elapsed;
    sf::Time diff;
    sf::Time period=sf::seconds(1.0/framerate);
    sf::Clock clock;

    std::vector<GameState*> states;
    states.push_back(new TitleScreen());
    sf::Vector2f rectsize;
    sf::Vector2f rectpos;
    sf::FloatRect bounds;
    bool ispressed=false;
    sf::Event event;
    sf::Vector2i mouse_pos;
    bool typing=false;
    bool hastyped=false;
    int typingindex=0;
    double dist=0.;
    std::string target;

    int controlindex=-1;
    bool wait=false;

    std::string name;
    std::string ip;
    unsigned short port;

    while (window.isOpen())
    {
        if (!wait)
        {
            mouse_pos=sf::Mouse::getPosition(window);
            mouse_pos.x=int(mouse_pos.x-0.5*sdiff);
            ispressed=sf::Mouse::isButtonPressed(sf::Mouse::Left);
            if (mouse_pos.x>=0 && mouse_pos.x<=dfactor*raw_width && mouse_pos.y>=0 && mouse_pos.y<=dfactor*raw_height)
            {
                for (int i=0;i<states.back()->_buttons.size();i++)
                {
                    rectpos=states.back()->_buttons[i]._shape.getPosition();
                    rectsize=states.back()->_buttons[i]._shape.getSize();
                    if (mouse_pos.x>=rectpos.x-0.5*rectsize.x && mouse_pos.x<=rectpos.x+0.5*rectsize.x
                        && mouse_pos.y>=rectpos.y-0.5*rectsize.y && mouse_pos.y<=rectpos.y+0.5*rectsize.y)
                    {
                        //touching the button.
                        if (ispressed)
                        {
                            states.back()->_buttons[i]._shape.setFillColor(states.back()->_buttons[i]._colour3);
                            states.back()->_buttons[i]._shape.setOutlineColor(states.back()->_buttons[i]._outlinecolour3);
                            states.back()->_buttons[i]._text.setFillColor(states.back()->_buttons[i]._textcolour3);
                        }
                        else if (states.back()->_buttons[i]._shape.getFillColor()==states.back()->_buttons[i]._colour3)
                        {
                            if (states.back()->_buttons[i]._controlchange==false)
                            {
                                //move onto the target state.
                                states.back()->_buttons[i]._shape.setFillColor(states.back()->_buttons[i]._colour1);
                                states.back()->_buttons[i]._shape.setOutlineColor(states.back()->_buttons[i]._outlinecolour1);
                                states.back()->_buttons[i]._text.setFillColor(states.back()->_buttons[i]._textcolour1);

                                target=states.back()->_buttons[i]._target;

                                if (target=="Quit") {delete states.back(); states.pop_back();}
                                else if (target=="Options") {states.push_back(new OptionsScreen());}
                                else if (target=="Singleplayer") {states.push_back(new SingleplayerScreen());}
                                else if (target=="Multiplayer") {states.push_back(new MultiplayerScreen(dfactor,std::string("192.168.823.106"),std::string("50000")));}
                                else if (target=="MultiplayerHost") {states.push_back(new GameScreen(dfactor,0,"",50000,name));}
                                else if (target=="MultiplayerJoin") {states.push_back(new GameScreen(dfactor,1,"",50000,name));}
                                else if (target=="SingleplayerAI") {states.push_back(new GameScreen(dfactor,2,"",50000,name));}
                                else if (target=="SingleplayerLineup") {states.push_back(new GameScreen(dfactor,3,"",50000,name));}
                                else if (target=="Controls") {states.push_back(new ControlScreen());}
                                else if (target=="Changecue") {states.push_back(new ChangeCueScreen());}
                                else if (target=="Default")
                                {
                                    user_controls=default_controls;
                                    sf::FloatRect bounds;
                                    for (int i=0;i<18;i++)
                                    {
                                        states.back()->_buttons[i]._text.setString(KeyToString(user_controls[states.back()->_buttons[i]._target]));
                                        bounds=states.back()->_buttons[i]._text.getLocalBounds();
                                        states.back()->_buttons[i]._text.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
                                    }
                                }
                                break;
                            }
                            else
                            {
                                controlindex=i;
                                states.back()->_buttons[controlindex]._text.setString("?");
                                bounds=states.back()->_buttons[controlindex]._text.getLocalBounds();
                                states.back()->_buttons[controlindex]._text.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
                                wait=true;
                                break;
                            }
                        }
                        else
                        {
                            states.back()->_buttons[i]._shape.setFillColor(states.back()->_buttons[i]._colour2);
                            states.back()->_buttons[i]._shape.setOutlineColor(states.back()->_buttons[i]._outlinecolour2);
                            states.back()->_buttons[i]._text.setFillColor(states.back()->_buttons[i]._textcolour2);
                        }
                    }
                    else
                    {
                        states.back()->_buttons[i]._shape.setFillColor(states.back()->_buttons[i]._colour1);
                        states.back()->_buttons[i]._shape.setOutlineColor(states.back()->_buttons[i]._outlinecolour1);
                        states.back()->_buttons[i]._text.setFillColor(states.back()->_buttons[i]._textcolour1);
                    }
                }
            }

            if (states.size()==0) {break;}

            if (mouse_pos.x>=0 && mouse_pos.x<=window_width && mouse_pos.y>=0 && mouse_pos.y<=window_height)
            {
                for (int i=0;i<states.back()->_inputboxes.size();i++)
                {
                    rectpos=states.back()->_inputboxes[i]._shape.getPosition();
                    rectsize=states.back()->_inputboxes[i]._shape.getSize();
                    if (ispressed)
                    {
                        if (mouse_pos.x>=rectpos.x-0.5*rectsize.x && mouse_pos.x<=rectpos.x+0.5*rectsize.x
                            && mouse_pos.y>=rectpos.y-0.5*rectsize.y && mouse_pos.y<=rectpos.y+0.5*rectsize.y)
                        {
                            states.back()->_inputboxes[i]._shape.setOutlineColor(states.back()->_inputboxes[i]._outlinecolour2);
                            sf::Color c=states.back()->_inputboxes[i]._backtext.getFillColor();
                            states.back()->_inputboxes[i]._backtext.setFillColor(sf::Color(c.r,c.g,c.b,0));
                            c=states.back()->_inputboxes[i]._cursor.getFillColor();
                            states.back()->_inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));

                            dist=999.;
                            for (int j=0;j<states.back()->_inputboxes[i]._input.length()+1;j++)
                            {
                                states.back()->_inputboxes[i]._text.setString(states.back()->_inputboxes[i]._input.substr(0,j));
                                bounds=states.back()->_inputboxes[i]._text.getGlobalBounds();
                                if (fabs(mouse_pos.x-bounds.left-bounds.width)<dist)
                                {
                                    dist=fabs(mouse_pos.x-bounds.left-bounds.width);
                                    states.back()->_inputboxes[i]._cursorpos=j;
                                    rectpos=states.back()->_inputboxes[i]._cursor.getPosition();
                                    states.back()->_inputboxes[i]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*states.back()->_inputboxes[typingindex]._abscursorthickness,rectpos.y));
                                }
                                states.back()->_inputboxes[i]._text.setString(states.back()->_inputboxes[i]._input);
                            }
                        }
                        else
                        {
                            states.back()->_inputboxes[i]._shape.setOutlineColor(states.back()->_inputboxes[i]._outlinecolour1);
                            sf::Color c=states.back()->_inputboxes[i]._backtext.getFillColor();
                            if (states.back()->_inputboxes[i]._input.empty())
                            {
                                states.back()->_inputboxes[i]._backtext.setFillColor(sf::Color(c.r,c.g,c.b,150));
                            }
                            c=states.back()->_inputboxes[i]._cursor.getFillColor();
                            states.back()->_inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,0));
                        }
                    }
                }
            }

            typing=false;
            for (int i=0;i<states.back()->_inputboxes.size();i++)
            {
                if (states.back()->_inputboxes[i]._shape.getOutlineColor()==states.back()->_inputboxes[i]._outlinecolour2)
                {
                    //typing in box. All inputs except mouse are ignored.
                    typing=true;
                    typingindex=i;
                    break;
                }
            }
        }

        hastyped=false;
        //deal with user input.
        while (window.pollEvent(event))
        {
            switch (event.type)
            {
                case sf::Event::Closed:
                    window.close();
                    break;
                case sf::Event::TextEntered:
                    if (!wait && typing && (event.text.unicode<127 && event.text.unicode>31))
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._input.insert(size_t(states.back()->_inputboxes[typingindex]._cursorpos),std::string(sf::String(event.text.unicode)));
                        states.back()->_inputboxes[typingindex]._cursorpos+=1;
                    }
                case sf::Event::KeyPressed:
                    if (!wait && typing && event.key.code==sf::Keyboard::Backspace)
                    {
                        hastyped=true;
                        if (states.back()->_inputboxes[typingindex]._cursorpos!=0)
                        {
                            states.back()->_inputboxes[typingindex]._input.erase(states.back()->_inputboxes[typingindex]._cursorpos-1,1);
                            states.back()->_inputboxes[typingindex]._cursorpos-=1;
                        }
                    }
                    else if (!wait && typing && event.key.code==sf::Keyboard::Left)
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._cursorpos=std::max(0,states.back()->_inputboxes[typingindex]._cursorpos-1);
                    }
                    else if (!wait && typing && event.key.code==sf::Keyboard::Right)
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._cursorpos=std::min(int(states.back()->_inputboxes[typingindex]._input.length()),states.back()->_inputboxes[typingindex]._cursorpos+1);
                    }
                    else if (!wait && typing && event.key.code==sf::Keyboard::Delete)
                    {
                        hastyped=true;
                        if (states.back()->_inputboxes[typingindex]._cursorpos!=states.back()->_inputboxes[typingindex]._input.length())
                        {
                            states.back()->_inputboxes[typingindex]._input.erase(states.back()->_inputboxes[typingindex]._cursorpos,1);
                        }
                    }
                    else if (wait && controlindex!=-1)
                    {
                        target=KeyToString(event.key.code);
                        user_controls[states.back()->_buttons[controlindex]._target]=event.key.code;
                        states.back()->_buttons[controlindex]._text.setString(target);
                        bounds=states.back()->_buttons[controlindex]._text.getLocalBounds();
                        states.back()->_buttons[controlindex]._text.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));

                        //reset the button colour.
                        states.back()->_buttons[controlindex]._shape.setFillColor(states.back()->_buttons[controlindex]._colour1);
                        states.back()->_buttons[controlindex]._shape.setOutlineColor(states.back()->_buttons[controlindex]._outlinecolour1);
                        states.back()->_buttons[controlindex]._text.setFillColor(states.back()->_buttons[controlindex]._textcolour1);

                        wait=false;
                        controlindex=-1;
                        break;
                    }
            }
        }

        if (!wait && hastyped)
        {
            states.back()->_inputboxes[typingindex]._t=0.;
            sf::Color c=states.back()->_inputboxes[typingindex]._cursor.getFillColor();
            states.back()->_inputboxes[typingindex]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));
            rectpos=states.back()->_inputboxes[typingindex]._cursor.getPosition();
            states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input.substr(0,states.back()->_inputboxes[typingindex]._cursorpos));
            bounds=states.back()->_inputboxes[typingindex]._text.getGlobalBounds();
            states.back()->_inputboxes[typingindex]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*states.back()->_inputboxes[typingindex]._abscursorthickness,rectpos.y));
            states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input);

            bounds=states.back()->_inputboxes[typingindex]._cursor.getGlobalBounds();
            rectpos=states.back()->_inputboxes[typingindex]._shape.getPosition();
            rectsize=states.back()->_inputboxes[typingindex]._shape.getSize();
            if (bounds.left+bounds.width>rectpos.x+0.5*rectsize.x)
            {
                states.back()->_inputboxes[typingindex]._input.resize(states.back()->_inputboxes[typingindex]._input.length()-1);
                if (states.back()->_inputboxes[typingindex]._cursorpos>states.back()->_inputboxes[typingindex]._input.length())
                {
                    states.back()->_inputboxes[typingindex]._cursorpos=states.back()->_inputboxes[typingindex]._input.length();
                }
                rectpos=states.back()->_inputboxes[typingindex]._cursor.getPosition();
                states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input.substr(0,states.back()->_inputboxes[typingindex]._cursorpos));
                bounds=states.back()->_inputboxes[typingindex]._text.getGlobalBounds();
                states.back()->_inputboxes[typingindex]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*states.back()->_inputboxes[typingindex]._abscursorthickness,rectpos.y));
                states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input);
            }
        }

        if (!wait) {states.back()->update(1.0/framerate);}

        window.clear(states.back()->_background);

        window.draw(p1);
        window.draw(p2);

        for (int i=0;i<states.back()->_shapes.size();i++)
        {
            window.draw(*states.back()->_shapes[i]);
        }
        for (int i=0;i<states.back()->_buttons.size();i++)
        {
            window.draw(states.back()->_buttons[i]._shape);
            window.draw(states.back()->_buttons[i]._text);
        }
        for (int i=0;i<states.back()->_inputboxes.size();i++)
        {
            window.draw(states.back()->_inputboxes[i]._shape);
            window.draw(states.back()->_inputboxes[i]._backtext);
            window.draw(states.back()->_inputboxes[i]._text);
            window.draw(states.back()->_inputboxes[i]._cursor);
        }
        window.display();

        elapsed=clock.restart();
        diff=period-elapsed;
        if (diff.asSeconds()>0) {sf::sleep(diff);}
    }
//    while (window.isOpen())
//    {
//        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space))
//        {
//            if (!placing_white && done)
//            {
//                //take the shot!
//                spin_selector_pos=spin_selector.getPosition();
//                spin_dot_pos=spin_dot.getPosition();
//                dist=(ball_radius/(0.4*panel_height))*sqrt(pow(spin_selector_pos.x-spin_dot_pos.x,2.)+pow(spin_selector_pos.y-spin_dot_pos.y,2.));
//                angle=atan2(spin_dot_pos.x-spin_selector_pos.x,spin_selector_pos.y-spin_dot_pos.y);
//                cue._offset=dist;
//                cue._theta=angle;
//                cue._speed=1.2*power;
//                //cue.perturb();
//                cue.shot();
//                balls[0]._vx=cue._ballv*sin(cue._angle);
//                balls[0]._vy=cue._ballv*cos(cue._angle);
//                balls[0]._xspin=cue._ballparspin*sin(cue._angle)+cue._ballperspin*cos(cue._angle);
//                balls[0]._yspin=cue._ballparspin*cos(cue._angle)-cue._ballperspin*sin(cue._angle);
//                balls[0]._rspin=cue._ballrspin;
//
//                packet.clear();
//                packetId=2;
//                packet << packetId << balls[0]._vx << balls[0]._vy << balls[0]._xspin << balls[0]._yspin << balls[0]._rspin;
//                socket.send(packet);
//                continue;
//            }
//        }
//        if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
//        {
//            if (!placing_white && done)
//            {
//                mouse_pos=sf::Mouse::getPosition(window);
//                if (mouse_pos.x>=0 && mouse_pos.x<=panel_height && mouse_pos.y>=window_height && mouse_pos.y<=window_height+panel_height)
//                {
//                    dist=fmin(sqrt(pow(mouse_pos.x-0.5*panel_height,2.)+pow(mouse_pos.y-0.5*panel_height-window_height,2.)),0.4*panel_height-panel_height*0.4*0.005/0.02625);
//                    angle=atan2(mouse_pos.x-0.5*panel_height,mouse_pos.y-0.5*panel_height-window_height);
//                    spin_dot.setPosition(sf::Vector2f(0.5*panel_height+dist*sin(angle),window_height+panel_height*0.5+dist*cos(angle)));
//
//                    spin_selector_pos=spin_selector.getPosition();
//                    spin_dot_pos=spin_dot.getPosition();
//                    dist=(spin_selector_pos.y-spin_dot_pos.y)*0.5;
//                    roots=qsolve_quadratic(1.,2*(sqrt(pow(0.2*panel_height,2.)-pow(dist,2.))*cos(alpha)+dist*sin(alpha)),pow(0.2*panel_height,2.)-pow(0.4*panel_height,2.));
//                    elevation_pointer.setPosition(sf::Vector2f(window_width-0.5*panel_height+sqrt(pow(0.2*panel_height,2.)-pow(dist,2.)),window_height+0.5*panel_height-dist));
//                    elevation_pointer.setRotation(-alpha*180./pi);
//                    if (roots[0]!=roots[0])
//                    {
//                        elevation_pointer.setSize(sf::Vector2f(roots[1],0.2*panel_height*0.005/0.02625));
//                    }
//                    else if (roots[1]!=roots[1])
//                    {
//                        elevation_pointer.setSize(sf::Vector2f(roots[0],0.2*panel_height*0.005/0.02625));
//                    }
//                    else
//                    {
//                        elevation_pointer.setSize(sf::Vector2f(fmax(roots[0],roots[1]),0.2*panel_height*0.005/0.02625));
//                    }
//                    change=true;
//                }
//            }
//        }
//
//        if (change)
//        {
//            spin_selector_pos=spin_selector.getPosition();
//            spin_dot_pos=spin_dot.getPosition();
//            dist=(ball_radius/(0.4*panel_height))*sqrt(pow(spin_selector_pos.x-spin_dot_pos.x,2.)+pow(spin_selector_pos.y-spin_dot_pos.y,2.));
//            angle=atan2(spin_dot_pos.x-spin_selector_pos.x,spin_selector_pos.y-spin_dot_pos.y);
//            cue._offset=dist;
//            cue._theta=angle;
//            cue._speed=1.2*power;
//            cue._alpha=alpha;
//        }

    return 0;
}
