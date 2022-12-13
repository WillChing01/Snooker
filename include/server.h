#ifndef SERVER_H_INCLUDED
#define SERVER_H_INCLUDED

class Server
{
    public:
        bool running=true;

        bool player_turn=0;

        int scores[2]={0,0};
        int frames[2]={0,0};
        int framesbestof=3;
        bool placing_white=true;
        bool touching=false;
        bool isfoul=false;
        bool ismiss=false;
        bool ispush=false;
        bool isredon=true;
        bool redsLeft=true;
        int colourClearOrder=2; //next colour when no reds left.
        int foulscore=0;
        int nom_colour_order=0; //2 yellow,3 green, etc.
        bool isfreeball=false;

        bool newIsRedOn=true;
        bool newIsFreeBall=false;

        std::array<double,66> posbefore;
        std::array<bool,22> potbefore;
        bool wasredon=false;
        bool wasfreeball=false;
        bool wasRedsLeft=false;

        bool gameover=false;
        int highbreak[2]={0,0};
        int centuries[2]={0,0};
        int current_break=0;

        std::vector<std::pair<std::string,std::string> > messageLog;

        //info for each shot.
        std::vector<int> ball_hit_order;
        std::vector<int> ball_potted_order;

        sf::IpAddress serverIp;
        unsigned short port;
        sf::TcpListener listener;
        bool _connectionError=false;

        std::array<sf::TcpSocket,2> players;
        std::array<sf::TcpSocket,4> spectators;

        std::array<std::string,2> pnames={"PLAYER 1","PLAYER 2"};
        std::array<std::string,4> snames;

        sf::Packet packet;
        sf::Uint16 packetId=0;

        //for objects just need balls and cushions.

        //prepare cushions.
        Cushion servercushions[6];

        Ball serverballs[22];

        //grid squares.
        std::array<std::array<std::array<int,7>,39>,73> grid={};
        std::array<std::array<int,39>,73> grid_index={};

        //matrix.
        Eigen::Matrix<double,46,46> sM_=Eigen::MatrixXd::Zero(46,46);

        std::vector<std::array<double,66> > result;

        Server(bool isOnline=true);
        std::vector<std::array<double,66> > simulate(Ball balls[22],Cushion cush[6]);
        Eigen::Matrix<double,46,1> collisions(Ball b[22], Cushion cush[6]);
        void handleIncomingConnections();
        void executionThread();
        void shutdown();
        void respot();
        void rackballs();
        void turnpacket();
        void sendNamePacket();
        void resetframe();
        bool is_ball_blocking(double x1, double y1, double x2, double y2, double x3, double y3);
        bool is_snookered(int ballOnOrder);
        void handleRedOn();
        void handleColourOn();
        void handleFreeBall();
        void applyRulesAndScore();
        void handleFoulMissResponse(std::string msg);
        void broadcast(std::string name, std::string message);
        void sendBallPositions();
        void sendShotSimulation();
        void resetServer();
};

Server::Server(bool isOnline)
{
    if (isOnline)
    {
        try {serverIp=sf::IpAddress::getLocalAddress();}
        catch (...) {_connectionError=true;}
        listener.setBlocking(false);
        if (listener.listen(sf::Socket::AnyPort)!=sf::Socket::Done)
        {
            //error.
            std::cout << "Tcp listener could not bind to port." << std::endl;
            _connectionError=true;
        }
        try{port=listener.getLocalPort();}
        catch (...) {_connectionError=true;}
    }

    for (int i=0;i<2;i++)
    {
        players[i].setBlocking(false);
    }
    for (int i=0;i<4;i++)
    {
        spectators[i].setBlocking(false);
    }

    servercushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
    servercushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
    servercushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
    servercushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
    servercushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
    servercushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

    //prepare matrix.
    for (int i=0;i<44;i++)
    {
        sM_(i,i)=1/ball_mass;
    }
    sM_(44,44)=0;
    sM_(45,45)=0;
}

void Server::broadcast(std::string name, std::string message)
{
    packet.clear();
    packetId=3;
    packet << packetId << name << message;
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            players[i].send(packet);
        }
    }
    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            spectators[i].send(packet);
        }
    }

    messageLog.push_back({name,message});
}

void Server::sendBallPositions()
{
    packet.clear();
    packet << sf::Uint16(1) << sf::Uint32(1);
    for (int i=0;i<22;i++)
    {
        if (serverballs[i]._potted)
        {
            packet << -100. << -100. << -100.;
        }
        else
        {
            packet << serverballs[i]._x;
            packet << serverballs[i]._y;
            packet << serverballs[i]._z;
        }
    }
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            players[i].send(packet);
        }
    }
    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            spectators[i].send(packet);
        }
    }
}

void Server::sendShotSimulation()
{
    //send the results in a packet.
    packet.clear();
    std::vector<double> flatresult;
    for (int i=0;i<result.size();i++)
    {
        for (int j=0;j<66;j++)
        {
            flatresult.push_back(result[i][j]);
        }
    }
    for (int i=0;i<22;i++)
    {
        if (!serverballs[i]._potted)
        {
            flatresult.push_back(serverballs[i]._x);
            flatresult.push_back(serverballs[i]._y);
            flatresult.push_back(serverballs[i]._z);
        }
        else
        {
            flatresult.push_back(-100.);
            flatresult.push_back(-100.);
            flatresult.push_back(-100.);
        }
    }
    packet << sf::Uint16(1) << sf::Uint32(result.size()+1);
    for (int i=0;i<flatresult.size();i++)
    {
        packet << flatresult[i];
    }

    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            players[i].send(packet);
        }
    }
    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            spectators[i].send(packet);
        }
    }
}

void Server::resetframe()
{
    scores[0]=0;
    scores[1]=0;

    if ((frames[0]+frames[1])%2==0) {player_turn=0;}
    else {player_turn=1;}

    isredon=true;
    isfreeball=false;
    placing_white=true;
    current_break=0;

    rackballs();

    sendBallPositions();
}

bool Server::is_ball_blocking(double x1, double y1, double x2, double y2, double x3, double y3)
{
    //1 and 2 are start/final positions of cueball.
    //3 is position of possible snookering ball.

    bool blocking=false;

    double dx1=x2-x1; double dy1=y2-y1;
    double dx2=x3-x1; double dy2=y3-y1;
    double dx3=x3-x2; double dy3=y3-y2;

    double a=std::acos((dx1*dx2+dy1*dy2)/(std::sqrt(dx1*dx1+dy1*dy1)*std::sqrt(dx2*dx2+dy2*dy2)));
    double b=std::acos((-dx1*dx3-dy1*dy3)/(std::sqrt(dx1*dx1+dy1*dy1)*std::sqrt(dx3*dx3+dy3*dy3)));

    if (a<=pi/2. && b<=pi/2.)
    {
        //possible snooker. check closest approach.
        double height=std::sqrt(dx2*dx2+dy2*dy2)*std::sin(a);
        if (height<2.*ball_radius) {blocking=true;}
        else {blocking=false;}
    }
    else
    {
        blocking=false;
    }

    return blocking;
}

bool Server::is_snookered(int ballOnOrder=8)
{
    bool snookered=false;
    bool leftHit=false;
    bool rightHit=false;
    double dx,dy;
    double px,py;
    double d;
    double x1,y1; double x2,y2;

    //check if the cue ball can hit both edges of the ball on.

    if (ballOnOrder>=8)
    {
        //red is on.
        leftHit=false;
        rightHit=false;
        for (int i=7;i<22;i++)
        {
            if (serverballs[i]._potted) {continue;}
            //check if snookered (reds cannot snooker each other).
            bool leftHit2=true;
            bool rightHit2=true;

            dx=serverballs[i]._x-serverballs[0]._x;
            dy=serverballs[i]._y-serverballs[0]._y;

            d=sqrt(dx*dx+dy*dy);
            px=-dy/d; py=dx/d;

            x1=serverballs[0]._x+dx*(1.-std::pow(2.*ball_radius/d,2.))+px*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));
            y1=serverballs[0]._y+dy*(1.-std::pow(2.*ball_radius/d,2.))+py*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));

            x2=serverballs[0]._x+dx*(1.-std::pow(2.*ball_radius/d,2.))-px*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));
            y2=serverballs[0]._y+dy*(1.-std::pow(2.*ball_radius/d,2.))-py*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));

            for (int j=1;j<7;j++)
            {
                if (serverballs[j]._potted) {continue;}

                //check if any side of red is accessible.
                bool isBlocking=is_ball_blocking(serverballs[0]._x,serverballs[0]._y,x1,y1,serverballs[j]._x,serverballs[j]._y);
                if (isBlocking) {leftHit2=false;}

                isBlocking=is_ball_blocking(serverballs[0]._x,serverballs[0]._y,x2,y2,serverballs[j]._x,serverballs[j]._y);
                if (isBlocking) {rightHit2=false;}
            }

            if (leftHit2) {leftHit=true;}
            if (rightHit2) {rightHit=true;}
        }

        if (leftHit && rightHit) {snookered=false;}
        else {snookered=true;}
    }
    else
    {
        //colour is on.
        dx=serverballs[ballOnOrder-1]._x-serverballs[0]._x;
        dy=serverballs[ballOnOrder-1]._y-serverballs[0]._y;

        d=sqrt(dx*dx+dy*dy);
        px=-dy/d; py=dx/d;

        x1=serverballs[0]._x+dx*(1.-std::pow(2.*ball_radius/d,2.))+px*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));
        y1=serverballs[0]._y+dy*(1.-std::pow(2.*ball_radius/d,2.))+py*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));

        x2=serverballs[0]._x+dx*(1.-std::pow(2.*ball_radius/d,2.))-px*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));
        y2=serverballs[0]._y+dy*(1.-std::pow(2.*ball_radius/d,2.))-py*2.*ball_radius*std::sqrt(1-std::pow(2.*ball_radius/d,2.));

        leftHit=true;
        rightHit=true;

        for (int i=1;i<22;i++)
        {
            if (serverballs[i]._order==ballOnOrder) {continue;}
            if (serverballs[i]._potted) {continue;}

            //check for possible snooker.
            bool isBlocking=is_ball_blocking(serverballs[0]._x,serverballs[0]._y,x1,y1,serverballs[i]._x,serverballs[i]._y);
            if (isBlocking) {leftHit=false;}

            isBlocking=is_ball_blocking(serverballs[0]._x,serverballs[0]._y,x2,y2,serverballs[i]._x,serverballs[i]._y);
            if (isBlocking) {rightHit=false;}
        }

        if (leftHit && rightHit) {snookered=false;}
        else {snookered=true;}
    }

    return snookered;
}

void Server::handleRedOn()
{
    //check that first ball hit is a red.
    if (ball_hit_order[0]<8) {isfoul=true; foulscore=std::max(foulscore,std::max(ball_hit_order[0],4));}

    //check which balls were potted.
    for (int i=0;i<ball_potted_order.size();i++)
    {
        if (ball_potted_order[i]<8) {isfoul=true; foulscore=std::max(foulscore,std::max(ball_potted_order[i],4));}
    }
    if (!isfoul)
    {
        //add points.
        scores[player_turn]+=ball_potted_order.size();
        current_break+=ball_potted_order.size();
        if (current_break>highbreak[player_turn]) {highbreak[player_turn]=current_break;}

        if (ball_potted_order.size()==0)
        {
            //no balls potted.
            current_break=0;
            player_turn=!player_turn;
            newIsRedOn=true;
            newIsFreeBall=false;
        }
        else
        {
            //potted a red.
            newIsRedOn=false;
            newIsFreeBall=false;
        }
    }
    else
    {
        //check if any reds left.
        if (redsLeft)
        {
            bool allPotted=true;
            for (int i=7;i<22;i++)
            {
                if (!serverballs[i]._potted) {allPotted=false; break;}
            }
            if (allPotted) {redsLeft=false;}

            if (!redsLeft)
            {
                newIsRedOn=false;
                newIsFreeBall=false;
            }
        }
    }
}

void Server::handleColourOn()
{
    //for colour on.
    //check that first ball hit is the correct colour.
    if (ball_hit_order[0]!=nom_colour_order)
    {
        isfoul=true;
        //ismiss???
        if (ball_hit_order[0]>=8)
        {
            foulscore=std::max(foulscore,4);
        }
        else
        {
            foulscore=std::max(foulscore,std::max(ball_hit_order[0],4));
        }
    }

    //check which balls were potted.
    bool colpot=false;
    for (int i=0;i<ball_potted_order.size();i++)
    {
        if (ball_potted_order[i]!=nom_colour_order)
        {
            isfoul=true;
            if (ball_potted_order[i]>=8)
            {
                foulscore=std::max(foulscore,4);
            }
            else
            {
                foulscore=std::max(foulscore,std::max(ball_potted_order[i],4));
            }
        }
        else
        {
            colpot=true;
        }
    }
    if (!isfoul)
    {
        //add points.
        if (colpot)
        {
            scores[player_turn]+=nom_colour_order;
            current_break+=nom_colour_order;
            if (current_break>highbreak[player_turn]) {highbreak[player_turn]=current_break;}

            if (!redsLeft) {colourClearOrder++;}
        }
        else
        {
            current_break=0;
            player_turn=!player_turn;
        }

        if (redsLeft)
        {
            newIsRedOn=true;
            newIsFreeBall=false;
        }
        else
        {
            newIsRedOn=false;
            newIsFreeBall=false;
        }
    }

    //check if any reds left.
    if (redsLeft)
    {
        bool allPotted=true;
        for (int i=7;i<22;i++)
        {
            if (!serverballs[i]._potted) {allPotted=false; break;}
        }
        if (allPotted) {redsLeft=false;}

        if (!redsLeft)
        {
            newIsRedOn=false;
            newIsFreeBall=false;
        }
    }
}

void Server::handleFreeBall()
{
    bool colhit=false;
    for (int i=0;i<ball_hit_order.size();i++)
    {
        if (ball_hit_order[i]<8)
        {
            if (ball_hit_order[i]!=nom_colour_order)
            {
                isfoul=true;
                foulscore=std::max(foulscore,std::max(ball_hit_order[i],4));
            }
            else {colhit=true;}
        }
    }
    if (!colhit) {isfoul=true; ismiss=true; foulscore=std::max(foulscore,4);}
    bool colpot=false;
    for (int i=0;i<ball_potted_order.size();i++)
    {
        if (ball_potted_order[i]>7) {colpot=true;}
        if (ball_potted_order[i]<8)
        {
            if (ball_potted_order[i]!=nom_colour_order)
            {
                isfoul=true;
                foulscore=std::max(foulscore,std::max(ball_potted_order[i],4));
            }
            else {colpot=true;}
        }
    }
    if (!isfoul)
    {
        //add points.
        if (colpot)
        {
            scores[player_turn]+=ball_potted_order.size();
            current_break+=ball_potted_order.size();
            if (current_break>highbreak[player_turn]) {highbreak[player_turn]=current_break;}
            newIsRedOn=false;
            newIsFreeBall=false;
        }
        else
        {
            current_break=0;
            player_turn=!player_turn;
            newIsRedOn=true;
            newIsFreeBall=false;

            //check if reds snookered by free ball.
            bool before=is_snookered();
            serverballs[nom_colour_order-1]._potted=true;
            bool after=is_snookered();
            serverballs[nom_colour_order-1]._potted=false;

            if (before && !after)
            {
                //snookered by the free ball.
            }
        }
    }
    else
    {
        scores[!player_turn]+=foulscore;
        current_break=0;
        player_turn=!player_turn;

        newIsFreeBall=is_snookered();
        if (newIsFreeBall) {newIsRedOn=false;}
        else {newIsRedOn=true;}
    }
}

void Server::applyRulesAndScore()
{
    isfoul=false;
    ismiss=false;
    ispush=false;
    foulscore=0;
    placing_white=false;
    //check if white potted.
    if (serverballs[0]._potted)
    {
        placing_white=true;
        isfoul=true;
        foulscore=std::max(foulscore,4);
        serverballs[0]._potted=false;
        serverballs[0]._x=cueball_break_x;
        serverballs[0]._y=cueball_break_y;
        serverballs[0]._z=ball_radius;
        serverballs[0]._vx=0.;
        serverballs[0]._vy=0.;
        serverballs[0]._vz=0.;
        serverballs[0]._xspin=0.;
        serverballs[0]._yspin=0.;
        serverballs[0]._rspin=0.;
    }

    //determine whether or not the shot was legal.

    newIsRedOn=false;
    newIsFreeBall=false;

    //was anything hit?
    if (ball_hit_order.size()==0)
    {
        isfoul=true; ismiss=true; foulscore=4;
        if (!isredon)
        {
            foulscore=std::max(nom_colour_order,4);
        }
    }
    else
    {
        if (isredon && redsLeft)
        {
            handleRedOn();
        }
        else if (!isfreeball || !redsLeft)
        {
            handleColourOn();
        }
        else if (isfreeball)
        {
            handleFreeBall();
        }
    }

    if (isfoul)
    {
        ismiss=true;
        respot(); // second time - if there was a foul on clearing the colours the colour is replaced.

        scores[!player_turn]+=foulscore;
        current_break=0;
        player_turn=!player_turn;

        if (redsLeft)
        {
            newIsFreeBall=is_snookered();
            if (newIsFreeBall) {newIsRedOn=false;}
            else {newIsRedOn=true;}
        }
        else
        {
            newIsFreeBall=is_snookered(nom_colour_order);
            newIsRedOn=false;
        }
    }

    isredon=newIsRedOn;
    isfreeball=newIsFreeBall;
}

void Server::handleFoulMissResponse(std::string msg)
{
    //replace balls etc.
    if (ismiss)
    {
        if (!msg.compare("/1"))
        {
            //play from here.
            isfoul=false;
            ismiss=false;
            turnpacket();
        }
        else if (!msg.compare("/2"))
        {
            //opponent play from here.
            isfoul=false;
            ismiss=false;
            isfreeball=false;
            player_turn=!player_turn;
            turnpacket();
        }
        else if (!msg.compare("/3"))
        {
            //replace balls for opponent.
            isfoul=false;
            ismiss=false;
            isfreeball=false;
            player_turn=!player_turn;

            for (int i=0;i<22;i++)
            {
                serverballs[i]._potted=potbefore[i];
                serverballs[i]._x=posbefore[3*i];
                serverballs[i]._y=posbefore[3*i+1];
                serverballs[i]._z=posbefore[3*i+2];
            }
            isredon=wasredon;
            isfreeball=wasfreeball;
            redsLeft=wasRedsLeft;

            sendBallPositions();

            turnpacket();
        }
    }
    else
    {
        //we know its a foul already.
        if (!msg.compare("/1"))
        {
            //play from here.
            isfoul=false;
            ismiss=false;
            turnpacket();
        }
        else if (!msg.compare("/2"))
        {
            //opponent play from here.
            isfoul=false;
            ismiss=false;
            player_turn=!player_turn;
            turnpacket();
        }
    }
}

void Server::turnpacket()
{
    //send turn packet.
    packetId=2;
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            packet << packetId;
            packet << (player_turn==i);
            packet << isfoul << ismiss << placing_white << isredon << isfreeball << gameover;
            packet << sf::Uint32(scores[i]) << sf::Uint32(scores[(i+1)%2]) << sf::Uint32(frames[i]) << sf::Uint32(frames[(i+1)%2]) << sf::Uint32(highbreak[i]) << sf::Uint32(highbreak[(i+1)%2]) << sf::Uint32(centuries[i]) << sf::Uint32(centuries[(i+1)%2]);
            packet << redsLeft << colourClearOrder;
            players[i].send(packet);
        }
    }
    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            packet << packetId;
            packet << false;
            packet << isfoul << ismiss << placing_white << isredon << isfreeball << gameover;
            packet << sf::Uint32(scores[0]) << sf::Uint32(scores[1]) << sf::Uint32(frames[0]) << sf::Uint32(frames[1]) << sf::Uint32(highbreak[0]) << sf::Uint32(highbreak[1]) << sf::Uint32(centuries[0]) << sf::Uint32(centuries[1]);
            packet << redsLeft << colourClearOrder;
            spectators[i].send(packet);
        }
    }
}

void Server::sendNamePacket()
{
    //send packet with names of playing players (not spectators).

    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            if (i==0)
            {
                //host.
                packet.clear();
                packetId=4;
                packet << packetId << pnames[0] << pnames[1] << framesbestof;
            }
            else
            {
                //joined player.
                packet.clear();
                packetId=4;
                packet << packetId << pnames[1] << pnames[0] << framesbestof;
            }
            players[i].send(packet);
        }
    }

    packet.clear();
    packetId=4;
    packet << packetId << pnames[0] << pnames[1] << framesbestof;

    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            spectators[i].send(packet);
        }
    }
}

void Server::rackballs()
{
    //cueball.
    serverballs[0]._x=cueball_break_x;
    serverballs[0]._y=cueball_break_y;
    serverballs[0]._order=1;
    //yellow.
    serverballs[1]._x=yellow_x;
    serverballs[1]._y=yellow_y;
    serverballs[1]._order=2;
    //green.
    serverballs[2]._x=green_x;
    serverballs[2]._y=green_y;
    serverballs[2]._order=3;
    //brown.
    serverballs[3]._x=brown_x;
    serverballs[3]._y=brown_y;
    serverballs[3]._order=4;
    //blue.
    serverballs[4]._x=blue_x;
    serverballs[4]._y=blue_y;
    serverballs[4]._order=5;
    //pink.
    serverballs[5]._x=pink_x;
    serverballs[5]._y=pink_y;
    serverballs[5]._order=6;
    //black.
    serverballs[6]._x=black_x;
    serverballs[6]._y=black_y;
    serverballs[6]._order=7;

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
            serverballs[i]._x=pink_x-0.1-2*ball_radius-row*sqrt(3.0)*(ball_radius+DOUBLE_EPSILON);
            serverballs[i]._y=pink_y+h*(ball_radius+DOUBLE_EPSILON);
            serverballs[i]._order=i+1;
            i+=1;
        }
    }

    for (int i=0;i<22;i++)
    {
        serverballs[i]._potted=false;
        serverballs[i]._z=ball_radius;
        serverballs[i]._vx=0.;
        serverballs[i]._vy=0.;
        serverballs[i]._vz=0.;
        serverballs[i]._xspin=0.;
        serverballs[i]._yspin=0.;
        serverballs[i]._rspin=0.;
        gxmin=int(floor((serverballs[i]._x-serverballs[i]._r)/(2.*serverballs[i]._r)));
        gxmax=int(floor((serverballs[i]._x+serverballs[i]._r)/(2.*serverballs[i]._r)));
        gymin=int(floor((serverballs[i]._y-serverballs[i]._r)/(2.*serverballs[i]._r)));
        gymax=int(floor((serverballs[i]._y+serverballs[i]._r)/(2.*serverballs[i]._r)));

        serverballs[i]._gpos[0][0]=gxmin;
        serverballs[i]._gpos[0][1]=gymin;
        serverballs[i]._gpos[1][0]=gxmax;
        serverballs[i]._gpos[1][1]=gymin;
        serverballs[i]._gpos[2][0]=gxmax;
        serverballs[i]._gpos[2][1]=gymax;
        serverballs[i]._gpos[3][0]=gxmin;
        serverballs[i]._gpos[3][1]=gymax;
    }
}

void Server::respot()
{
    double xspot=0.;
    double yspot=0.;
    std::array<bool,6> covered;
    double dx=0.001/in2m;

    for (int i=0;i<6;i++)
    {
        covered[i]=false;
        xspot=colourpos[i][0];
        yspot=colourpos[i][1];
        serverballs[i+1]._z=ball_radius;

        for (int j=0;j<22;j++)
        {
            if (sqrt(pow(xspot-serverballs[j]._x,2.)+pow(yspot-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
            {
                covered[i]=true;
                break;
            }
        }
    }

    for (int i=6;i>0;i--)
    {
        if (!serverballs[i]._potted) {continue;}
        //do not respot if clearing the colours and not a foul.
        if (!redsLeft && serverballs[i]._order==colourClearOrder && !isfoul) {continue;}
        if (!covered[i-1])
        {
            serverballs[i]._x=colourpos[i-1][0];
            serverballs[i]._y=colourpos[i-1][1];
            covered[i-1]=true;
            serverballs[i]._potted=false;
            continue;
        }
        //default spot covered.
        //find highest empty spot.
        for (int j=5;j>-1;j--)
        {
            if (!covered[j])
            {
                serverballs[i]._x=colourpos[j][0];
                serverballs[i]._y=colourpos[j][1];
                covered[j]=true;
                serverballs[i]._potted=false;
                break;
            }
        }
        if (!serverballs[i]._potted) {continue;}

        //move closest to its spot.
        serverballs[i]._x=colourpos[i-1][0];
        serverballs[i]._y=colourpos[i-1][1];
        bool good=true;
        while (true)
        {
            serverballs[i]._x-=dx;
            if (serverballs[i]._x-rail_thickness-cush_thickness<epsilon)
            {
                serverballs[i]._x=-100.;
                break;
            }
            good=true;
            for (int j=0;j<22;j++)
            {
                if (i==j) {continue;}
                if (sqrt(pow(serverballs[i]._x-serverballs[j]._x,2.)+pow(serverballs[i]._y-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
                {
                    good=false;
                    break;
                }
            }
            if (good) {break;}
        }
        if (!serverballs[i]._potted) {continue;}

        serverballs[i]._x=colourpos[i-1][0];
        while (true)
        {
            serverballs[i]._x+=dx;
            good=true;
            for (int j=0;j<22;j++)
            {
                if (i==j) {continue;}
                if (sqrt(pow(serverballs[i]._x-serverballs[j]._x,2.)+pow(serverballs[i]._y-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
                {
                    good=false;
                    break;
                }
            }
            if (good) {break;}
        }
    }
}

Eigen::Matrix<double,46,1> Server::collisions(Ball b[22],Cushion cush[6])
{
    //broad phase collision check.
    Eigen::MatrixXd dv=Eigen::MatrixXd::Zero(46,1);

    int val;
    int thing;
    double dx;
    double dy;
    double dist;
    double relspeed;

    std::set<std::tuple<int,int>> collided;
    std::vector<int> already;

    Eigen::MatrixXd Z(46,0);
    int col=0;

    for (int i=0;i<22;i++)
    {
        if (b[i]._potted==false)
        {
            if (b[i]._vx!=0. || b[i]._vy!=0. || b[i]._vz!=0.)
            {
                already.clear();
                for (int j=0;j<4;j++)
                {
                    val=grid_index[b[i]._gpos[j][0]][b[i]._gpos[j][1]];
                    if (val>1)
                    {
                        for (int k=0;k<val;k++)
                        {
                            //add to potentially touching list.
                            thing=grid[b[i]._gpos[j][0]][b[i]._gpos[j][1]][k];
                            if (thing!=b[i]._order && b[thing-1]._potted==false)
                            {
                                if (std::find(already.begin(),already.end(),thing)==already.end())
                                {
                                    if (i>(thing-1) && collided.find(std::make_tuple(thing-1,i))!=collided.end())
                                    {
                                        continue;
                                    }
                                    else if (i<(thing-1) && collided.find(std::make_tuple(i,thing-1))!=collided.end())
                                    {
                                        continue;
                                    }
                                    already.push_back(thing);
                                    //get the distance using pythag.
                                    dx=b[thing-1]._x-b[i]._x;
                                    dy=b[thing-1]._y-b[i]._y;
                                    dist=sqrt(pow(dx,2.)+pow(dy,2.));
                                    relspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(dx,dy)-atan2(b[i]._vx,b[i]._vy))-sqrt(pow(b[thing-1]._vx,2.)+pow(b[thing-1]._vy,2.))*cos(atan2(dx,dy)-atan2(b[thing-1]._vx,b[thing-1]._vy));

                                    if (dist-2.*ball_radius<epsilon && relspeed>0.)
                                    {
                                        //actual collision.
                                        if (i>(thing-1))
                                        {
                                            collided.insert(std::make_tuple(thing-1,i));
                                        }
                                        else
                                        {
                                            collided.insert(std::make_tuple(i,thing-1));
                                        }
                                        Z.conservativeResizeLike(Eigen::MatrixXd::Zero(46,col+1));
                                        Z(2*i,col)=-dx/dist;
                                        Z(2*i+1,col)=-dy/dist;
                                        Z(2*(thing-1),col)=dx/dist;
                                        Z(2*(thing-1)+1,col)=dy/dist;
                                        col+=1;
                                    }
                                }
                            }
                        }
                    }
                }
                //check cushion collisions.
                if (fabs(blue_x-b[i]._x)>table_width-ball_radius-2.*epsilon || fabs(blue_y-b[i]._y)>0.5*table_width-ball_radius-2.*epsilon)
                {
                    double min_dist=999.;
                    double normal_angle;
                    double angle;
                    for (int j=0;j<6;j++)
                    {
                        std::tie(dist,angle)=cush[j].distance(b[i]._x,b[i]._y,b[i]._z);
                        if (dist<min_dist)
                        {
                            min_dist=dist;
                            normal_angle=angle;
                        }
                    }
                    relspeed=-sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(normal_angle-atan2(b[i]._vx,b[i]._vy));
                    if (min_dist-ball_radius<epsilon && relspeed>0.)
                    {
                        //collision with cushion.
                        Z.conservativeResizeLike(Eigen::MatrixXd::Zero(46,col+1));
                        Z(2*i,col)=sin(normal_angle);
                        Z(2*i+1,col)=cos(normal_angle);
                        Z(44,col)=-sin(normal_angle);
                        Z(45,col)=-cos(normal_angle);
                        col+=1;
                        //adjust the vertical spin off cushion.
                        double relspin=b[i]._xspin*sin(normal_angle+pi)+b[i]._yspin*cos(normal_angle+pi);
                        double dw;
                        if (relspin>0.)
                        {
                            //topspin.
                            dw=-5.*relspeed*(muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin<0.)
                        {
                            //backspin.
                            dw=-5.*relspeed*(-muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin==0.)
                        {
                            dw=-5.*relspeed*cushion_diff/pow(ball_radius,2.);
                        }
                        dw=fmax(dw,-relspin-relspeed/ball_radius);
                        b[i]._xspin+=dw*sin(normal_angle+pi);
                        b[i]._yspin+=dw*cos(normal_angle+pi);
                        //adjust sidespin off cushion.
                        double parspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*sin(atan2(b[i]._vx,b[i]._vy)-normal_angle);
                        double dvpar;

                        if ((parspeed>=0. && b[i]._rspin*ball_radius>parspeed) || (parspeed<0. && b[i]._rspin*ball_radius>parspeed))
                        {
                            dw=-5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius<parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        else if ((parspeed>=0. && b[i]._rspin*ball_radius<parspeed) || (parspeed<0. && b[i]._rspin*ball_radius<parspeed))
                        {
                            dw=5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius>parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        b[i]._rspin+=dw;
                        double _speed;
                        double _angle;
                        std::tie(_speed,_angle)=add_vectors(sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.)),atan2(b[i]._vx,b[i]._vy),dvpar,normal_angle+0.5*pi);
                        b[i]._vx=_speed*sin(_angle);
                        b[i]._vy=_speed*cos(_angle);
                    }
                }
            }
        }
    }

    //simultaneous.
    if (col>0)
    {
        Eigen::MatrixXd v(46,1);
        Eigen::MatrixXd J(col,1);

        for (int i=0;i<22;i++)
        {
            v(2*i,0)=b[i]._vx;
            v(2*i+1,0)=b[i]._vy;
        }
        v(44,0)=0;
        v(45,0)=0;

        //J=(Z.transpose()*M_*Z).colPivHouseholderQr().solve(-(1.+eball)*Z.transpose()*v);
        J=(Z.transpose()*sM_*Z).householderQr().solve(-(1.+eball)*Z.transpose()*v);
        dv=sM_*Z*J;
    }
    return dv;
}

std::vector<std::array<double,66> > Server::simulate(Ball balls[22],Cushion cush[6])
{
    //assumes that the input is not a static scenario (where all balls are still).
    std::vector<std::array<double,66> > pos;
    std::array<double,66> temp;
    std::array<double,3> out;
    std::array<double,3> out2;
    std::array<double,3> out3;
    std::array<double,3> out4;
    std::array<double,4> quartic;
    std::array<double,2> quadratic;
    std::array<std::array<double,5>,6> ctemp;
    std::array<double,5> x2;
    std::array<double,5> y2;
    std::array<double,5> z2;
    std::array<double,9> toctic;
    std::array<double,8> monoctic;
    std::array<double,8> octic;
    Eigen::Matrix<double,46,1> dv;
    Eigen::Matrix<double,46,1> zero;

    for (int i=0;i<46;i++)
    {
        dv(i,0)=0.;
        zero(i,0)=0.;
    }

    //set up the equations for each ball.

    double t; //time until next event.
    double totaltime=0.; //total time elapsed.
    double start;
    int c;

    //for solving quartics.
    double a1;
    double b1;
    double c1;
    double a2;
    double b2;
    double c2;
    double a3;
    double b3;
    double c3;
    double t0;
    double t1;
    double t2;
    double t3;
    double t4;
    double x;
    double y;
    double phi;
    double nspeed;
    double parspeed;
    double vtheta;
    double vphi;
    double normal;

    double xmin;
    double xmax;
    double xmin2;
    double xmax2;
    double ymin;
    double ymax;
    double ymin2;
    double ymax2;

    int xpos;
    int ypos;
    int gxmin;
    int gxmax;
    int gymin;
    int gymax;

    bool collision=false;

    ball_hit_order.clear();
    ball_potted_order.clear();
    bool hitdone=false;

    while (true)
    {
//        std::cout << "Time: " << totaltime << std::endl;

//        if (totaltime>0.2)
//        {
//            break;
//        }

        collision=false;
        t=999.;
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted==false)
            {
                balls[i].update_equation();
                t=fmin(t,balls[i]._t);
            }
        }
        if (t>998.)
        {
            //add final positions.
            for (int i=0;i<22;i++)
            {
                if (!balls[i]._potted)
                {
                    temp[3*i]=balls[i]._x;
                    temp[3*i+1]=balls[i]._y;
                    temp[3*i+2]=balls[i]._z;
                }
                else
                {
                    temp[3*i]=-100.;
                    temp[3*i+1]=-100.;
                    temp[3*i+2]=-100.;
                }
            }
            pos.push_back(temp);
            break;
        }
        //check for collisions (cushion/balls) or if ball into pocket.
        //check moving balls only.
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted || balls[i]._onrim)
            {
                continue;
            }
            if (balls[i]._vx==0. && balls[i]._vy==0. && balls[i]._vz==0. && balls[i]._xspin==0. && balls[i]._yspin==0.)
            {
                continue;
            }

            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmin-=1.01*ball_radius;
            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmax+=1.01*ball_radius;
            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymin-=1.01*ball_radius;
            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymax+=1.01*ball_radius;

            //check for ball-ball collisions.
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted || balls[j]._onrim || i==j)
                {
                    continue;
                }

                //initial check for i.
                xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmin2-=1.01*ball_radius;
                xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmax2+=1.01*ball_radius;
                ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymin2-=1.01*ball_radius;
                ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymax2+=1.01*ball_radius;

                if (ymin>ymax2 || ymin2>ymax || xmin>xmax2 || xmin2>xmax)
                {
                    continue;
                }
                //check if their trajectories are ever 2r apart.
                a1=balls[i]._ax[2]-balls[j]._ax[2];
                b1=balls[i]._ax[1]-balls[j]._ax[1];
                c1=balls[i]._ax[0]-balls[j]._ax[0];
                a2=balls[i]._ay[2]-balls[j]._ay[2];
                b2=balls[i]._ay[1]-balls[j]._ay[1];
                c2=balls[i]._ay[0]-balls[j]._ay[0];
                a3=balls[i]._az[2]-balls[j]._az[2];
                b3=balls[i]._az[1]-balls[j]._az[1];
                c3=balls[i]._az[0]-balls[j]._az[0];

                t0=pow(c1,2.)+pow(c2,2.)+pow(c3,2.)-4.*pow(ball_radius,2.);
                t1=2.*(b1*c1+b2*c2+b3*c3);
                t2=2.*(a1*c1+a2*c2+a3*c3)+pow(b1,2.)+pow(b2,2.)+pow(b3,2.);
                t3=2.*(a1*b1+a2*b2+a3*b3);
                t4=pow(a1,2.)+pow(a2,2.)+pow(a3,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        //verify the time of actual collision to be certain.
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                        out2[0]=balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0];
                        out2[1]=balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0];
                        out2[2]=balls[j]._az[2]*t0*t0+balls[j]._az[1]*t0+balls[j]._az[0];
                        out3[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        out3[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        out3[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
                        out4[0]=balls[j]._avx[1]*t0+balls[j]._avx[0];
                        out4[1]=balls[j]._avy[1]*t0+balls[j]._avy[0];
                        out4[2]=balls[j]._avz[1]*t0+balls[j]._avz[0];
                        a1=get_relspeed(out,out3,out2,out4);
                        if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                        {
                            //collision at the specified time.
                            t=t0;
                            collision=true;

                            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmin-=1.01*ball_radius;
                            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmax+=1.01*ball_radius;
                            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymin-=1.01*ball_radius;
                            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymax+=1.01*ball_radius;
                        }
                        else
                        {
                            //adjust the time minutely to ensure a collision.
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ax[1]*(t0-c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ay[1]*(t0-c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._az[1]*(t0-c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0-c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0-c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0-c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ax[1]*(t0+c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ay[1]*(t0+c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._az[1]*(t0+c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0+c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0+c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0+c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if(t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

//            //check if ball hits the baize.
//            if (balls[i]._inflight)
//            {
//                quadratic=qsolve_quadratic(balls[i]._az[2],balls[i]._az[1],balls[i]._az[0]-ball_radius);
//                for (int j=0;j<2;j++)
//                {
//                    if (quadratic[j]==quadratic[j] && quadratic[j]>=0. && quadratic[j]<t)
//                    {
//                        t0=quadratic[j];
//                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                        out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                        out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                        //check if not in a pocket.
//                        a1=999.;
//                        a2=999.;
//                        for (int k=0;k<4;k++)
//                        {
//                            //corner check
//                            a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                        }
//                        for (int k=0;k<2;k++)
//                        {
//                            a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                        }
//                        if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                        {
//                            t=t0;
//                            collision=true;
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                //check if not in a pocket.
//                                a1=999.;
//                                a2=999.;
//                                for (int k=0;k<4;k++)
//                                {
//                                    //corner check
//                                    a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                                }
//                                for (int k=0;k<2;k++)
//                                {
//                                    a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                                }
//                                if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                //check if not in a pocket.
//                                a1=999.;
//                                a2=999.;
//                                for (int k=0;k<4;k++)
//                                {
//                                    //corner check
//                                    a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                                }
//                                for (int k=0;k<2;k++)
//                                {
//                                    a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                                }
//                                if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }

            //broad check for cushions - eliminate doomed cases.
            if (xmin>xlim[0] && xmax<xlim[1] && ymin>ylim[0] && ymax<ylim[1])
            {
                continue;
            }

            //check for collision with cushions.

            //straight x.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ax[2],balls[i]._ax[1],balls[i]._ax[0]-xlim[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //straight y.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ay[2],balls[i]._ay[1],balls[i]._ay[0]-ylim[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        if (j==0) {a2=-out2[1];}
                        else {a2=out2[1];}
                        if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //check if hits rim.
            if (balls[i]._inflight)
            {
//                for (int j=0;j<2;j++)
//                {
//                    //solve the octic(!)
//                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
//                    x2[3]=2.*balls[i]._ax[2]*balls[i]._ax[1];
//                    x2[2]=2.*balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
//                    x2[1]=2.*balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0]);
//                    x2[0]=(balls[i]._ax[0]-mpockets[j][0])*(balls[i]._ax[0]-mpockets[j][0]);
//                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
//                    y2[3]=2.*balls[i]._ay[2]*balls[i]._ay[1];
//                    y2[2]=2.*balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
//                    y2[1]=2.*balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]);
//                    y2[0]=(balls[i]._ay[0]-mpockets[j][1])*(balls[i]._ay[0]-mpockets[j][1]);
//                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
//                    z2[3]=2.*balls[i]._az[2]*balls[i]._az[1];
//                    z2[2]=2.*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
//                    z2[1]=2.*balls[i]._az[1]*balls[i]._az[0];
//                    z2[0]=balls[i]._az[0]*balls[i]._az[0];
//
//                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
//                    toctic[7]=(2.*x2[4]*x2[3])+(2.*y2[4]*y2[3])+(2.*z2[4]*z2[3]);
//                    toctic[6]=(2.*x2[4]*x2[2]+x2[3]*x2[3])+(2.*y2[4]*y2[2]+y2[3]*y2[3])+(2.*z2[4]*z2[2]+z2[3]*z2[3]);
//                    toctic[5]=(2.*x2[4]*x2[1]+2.*x2[3]*x2[2])+(2.*y2[4]*y2[1]+2*y2[3]*y2[2])+(2.*z2[4]*z2[1]+2.*z2[3]*z2[2]);
//                    toctic[4]=(2.*x2[4]*x2[0]+2.*x2[3]*x2[1]+x2[2]*x2[2])+(2.*y2[4]*y2[0]+2.*y2[3]*y2[1]+y2[2]*y2[2])+(2.*z2[4]*z2[0]+2.*z2[3]*z2[1]+z2[2]*z2[2]);
//                    toctic[3]=(2.*x2[3]*x2[0]+2.*x2[2]*x2[1])+(2.*y2[3]*y2[0]+2.*y2[2]*y2[1])+(2.*z2[3]*z2[0]+2.*z2[2]*z2[1]);
//                    toctic[2]=(2.*x2[2]*x2[0]+x2[1]*x2[1])+(2.*y2[2]*y2[0]+y2[1]*y2[1])+(2.*z2[2]*z2[0]+z2[1]*z2[1]);
//                    toctic[1]=(2.*x2[1]*x2[0])+(2.*y2[1]*y2[0])+(2.*z2[1]*z2[0]);
//                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);
//
//                    toctic[8]+=2.*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
//                    toctic[7]+=2.*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
//                    toctic[6]+=2.*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
//                    toctic[5]+=2.*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
//                    toctic[4]+=2.*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
//                    toctic[3]+=2.*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
//                    toctic[2]+=2.*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
//                    toctic[1]+=2.*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
//                    toctic[0]+=2.*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));
//
//                    toctic[4]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[4]+y2[4])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[4];
//                    toctic[3]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[3]+y2[3])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[3];
//                    toctic[2]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[2]+y2[2])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[2];
//                    toctic[1]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[1]+y2[1])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[1];
//                    toctic[0]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[0]+y2[0])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[0];
//
//                    toctic[0]+=pow(pow(mpocket_r,2.)-pow(ball_radius,2.),2.);
//
//                    for (int k=0;k<8;k++)
//                    {
//                        monoctic[k]=toctic[k]/toctic[8];
//                    }
//                    octic=qsolve_octic(monoctic);
//                    //now to check the roots.
//                    for (int k=0;k<8;k++)
//                    {
//                        if (octic[k]==octic[k] && octic[k]>=0. && octic[k]<t)
//                        {
//                            t0=octic[k];
//                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
//                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
//                            out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                            a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                            phi=acos(out[2]/a1);
//                            nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                            vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                            if (a1<ball_radius && normal>0.)
//                            {
//                                t=t0;
//                                collision=true;
//                            }
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//                for (int j=0;j<4;j++)
//                {
//                    //solve the octic(!)
//                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
//                    x2[3]=2.*balls[i]._ax[2]*balls[i]._ax[1];
//                    x2[2]=2.*balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
//                    x2[1]=2.*balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0]);
//                    x2[0]=(balls[i]._ax[0]-cpockets[j][0])*(balls[i]._ax[0]-cpockets[j][0]);
//                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
//                    y2[3]=2.*balls[i]._ay[2]*balls[i]._ay[1];
//                    y2[2]=2.*balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
//                    y2[1]=2.*balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]);
//                    y2[0]=(balls[i]._ay[0]-cpockets[j][1])*(balls[i]._ay[0]-cpockets[j][1]);
//                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
//                    z2[3]=2.*balls[i]._az[2]*balls[i]._az[1];
//                    z2[2]=2.*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
//                    z2[1]=2.*balls[i]._az[1]*balls[i]._az[0];
//                    z2[0]=balls[i]._az[0]*balls[i]._az[0];
//
//                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
//                    toctic[7]=(2.*x2[4]*x2[3])+(2.*y2[4]*y2[3])+(2.*z2[4]*z2[3]);
//                    toctic[6]=(2.*x2[4]*x2[2]+x2[3]*x2[3])+(2.*y2[4]*y2[2]+y2[3]*y2[3])+(2.*z2[4]*z2[2]+z2[3]*z2[3]);
//                    toctic[5]=(2.*x2[4]*x2[1]+2.*x2[3]*x2[2])+(2.*y2[4]*y2[1]+2.*y2[3]*y2[2])+(2.*z2[4]*z2[1]+2.*z2[3]*z2[2]);
//                    toctic[4]=(2.*x2[4]*x2[0]+2.*x2[3]*x2[1]+x2[2]*x2[2])+(2.*y2[4]*y2[0]+2.*y2[3]*y2[1]+y2[2]*y2[2])+(2.*z2[4]*z2[0]+2.*z2[3]*z2[1]+z2[2]*z2[2]);
//                    toctic[3]=(2.*x2[3]*x2[0]+2.*x2[2]*x2[1])+(2.*y2[3]*y2[0]+2.*y2[2]*y2[1])+(2.*z2[3]*z2[0]+2.*z2[2]*z2[1]);
//                    toctic[2]=(2.*x2[2]*x2[0]+x2[1]*x2[1])+(2.*y2[2]*y2[0]+y2[1]*y2[1])+(2.*z2[2]*z2[0]+z2[1]*z2[1]);
//                    toctic[1]=(2.*x2[1]*x2[0])+(2.*y2[1]*y2[0])+(2.*z2[1]*z2[0]);
//                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);
//
//                    toctic[8]+=2.*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
//                    toctic[7]+=2.*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
//                    toctic[6]+=2.*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
//                    toctic[5]+=2.*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
//                    toctic[4]+=2.*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
//                    toctic[3]+=2.*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
//                    toctic[2]+=2.*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
//                    toctic[1]+=2.*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
//                    toctic[0]+=2.*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));
//
//                    toctic[4]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[4]+y2[4])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[4];
//                    toctic[3]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[3]+y2[3])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[3];
//                    toctic[2]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[2]+y2[2])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[2];
//                    toctic[1]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[1]+y2[1])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[1];
//                    toctic[0]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[0]+y2[0])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[0];
//
//                    toctic[0]+=pow(pow(cpocket_r,2.)-pow(ball_radius,2.),2.);
//
//                    for (int k=0;k<8;k++)
//                    {
//                        monoctic[k]=toctic[k]/toctic[8];
//                    }
//                    octic=qsolve_octic(monoctic);
//                    //now to check the roots.
//                    for (int k=0;k<8;k++)
//                    {
//                        if (octic[k]==octic[k] && octic[k]>=0. && octic[k]<t)
//                        {
//                            t0=octic[k];
//                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
//                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
//                            out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                            a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                            phi=acos(out[2]/a1);
//                            nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                            vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                            if (a1<ball_radius && normal>0.)
//                            {
//                                t=t0;
//                                collision=true;
//                            }
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
            }
            else
            {
                for (int j=0;j<2;j++)
                {
                    //middle pockets.
                    t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-mpockets[j][0],2.)+pow(balls[i]._ay[0]-mpockets[j][1],2.)-pow(mpocket_r,2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                            a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
                            if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                t=t0;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                                a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                                a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
                for (int j=0;j<4;j++)
                {
                    //corner pockets.
                    t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-cpockets[j][0],2.)+pow(balls[i]._ay[0]-cpockets[j][1],2.)-pow(cpocket_r,2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                            a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
                            if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                t=t0;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                                a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                                a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //check curved pocket cushions.
            for (int j=0;j<6;j++)
            {
                if (cush[j]._p1==0)
                {
                    ctemp[0]=cp[0];
                    ctemp[1]=cp[1];
                    ctemp[2]=cp[2];
                }
                else
                {
                    ctemp[0]=mp[0];
                    ctemp[1]=mp[1];
                    ctemp[2]=mp[2];
                }
                if (cush[j]._p2==0)
                {
                    ctemp[3]=cp[0];
                    ctemp[4]=cp[1];
                    ctemp[5]=cp[2];
                }
                else
                {
                    ctemp[3]=mp[0];
                    ctemp[4]=mp[1];
                    ctemp[5]=mp[2];
                }
                for (int k=3;k<6;k++)
                {
                    ctemp[k][0]=cush[j]._length-ctemp[k][0];
                    t0=ctemp[k][3];
                    ctemp[k][3]=-ctemp[k][4];
                    ctemp[k][4]=-t0;
                }
                for (int k=0;k<6;k++)
                {
                    x=cush[j]._x+(ctemp[k][0]*cos(cush[j]._angle)+ctemp[k][1]*sin(cush[j]._angle));
                    y=cush[j]._y+(-ctemp[k][0]*sin(cush[j]._angle)+ctemp[k][1]*cos(cush[j]._angle));
                    t4=pow(balls[i]._ax[2],2.)+pow(balls[i]._ay[2],2.);
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-x)+balls[i]._ay[2]*(balls[i]._ay[0]-y))+pow(balls[i]._ax[1],2.)+pow(balls[i]._ay[1],2.);
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-x)+balls[i]._ay[1]*(balls[i]._ay[0]-y));
                    t0=pow(balls[i]._ax[0]-x,2.)+pow(balls[i]._ay[0]-y,2.)-pow(ctemp[k][2],2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    //check the validity of quartic solutions.
                    for (int l=0;l<4;l++)
                    {
                        if (quartic[l]==quartic[l] && quartic[l]>=0. && quartic[l]<t)
                        {
                            //verify the time of collision.
                            t0=quartic[l];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                            a1=atan2(out[0]-x,out[1]-y);
                            std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                            a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                            a3=ctemp[k][3]+cush[j]._angle;
                            a3=a3-2.*pi*floor(a3/(2.*pi));
                            b3=ctemp[k][4]+cush[j]._angle;
                            b3=b3-2.*pi*floor(b3/(2.*pi));
                            if (a3>pi) {a3=a3-2.*pi;}
                            if (a3<-pi) {a3=a3+2.*pi;}
                            if (b3>pi) {b3=b3-2.*pi;}
                            if (b3<-pi) {b3=b3+2.*pi;}
                            if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
                            {
                                t=t0;
                                collision=true;
                            }
                            else
                            {
                                c=0;
                                while (true)
                                {
                                    c+=1;
                                    out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
                                    {
                                        if (t0-c*epsilon>=0.)
                                        {
                                            t=t0-c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
                                    {
                                        if (t0+c*epsilon<t)
                                        {
                                            t=t0+c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    if (c>1000)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //check middle straight bits.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ax[2],balls[i]._ax[1],balls[i]._ax[0]-mpline[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //check corner straight bits.
            for (int j=0;j<8;j++)
            {
                a1=0.5*balls[i]._ax[2]-0.5*cpline[j][0]*balls[i]._ay[2];
                a2=0.5*balls[i]._ax[1]-0.5*cpline[j][0]*balls[i]._ay[1];
                a3=0.5*balls[i]._ax[0]-0.5*cpline[j][0]*balls[i]._ay[0]+0.5*cpline[j][0]*cpline[j][1];
                b1=0.5*balls[i]._ay[2]-0.5*cpline[j][0]*balls[i]._ax[2];
                b2=0.5*balls[i]._ay[1]-0.5*cpline[j][0]*balls[i]._ax[1];
                b3=0.5*balls[i]._ay[0]-0.5*cpline[j][0]*balls[i]._ax[0]-0.5*cpline[j][1];

                t4=pow(a1,2.)+pow(b1,2.);
                t3=2.*(a1*a2+b1*b2);
                t2=2.*(a1*a3+b1*b3)+pow(a2,2.)+pow(b2,2.);
                t1=2.*(a2*a3+b2*b3);
                t0=pow(a3,2.)+pow(b3,2.)-pow(ball_radius,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                        b1=cpline[j][0]*a1+cpline[j][1];
                        a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                        a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                        if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i=0;i<22;i++)
        {
            //deal with balls on the rim.
            if (!balls[i]._onrim || balls[i]._potted)
            {
                continue;
            }

            //check when ball hits cushion rim.
            //work out which pocket ball is rolling in.
            c=0;
            for (int j=0;j<2;j++)
            {
                a1=sqrt(pow(balls[i]._x-mpockets[j][0],2.)+pow(balls[i]._y-mpockets[j][1],2.));
                if (a1<mpocket_r+epsilon)
                {
                    //rolling in this pocket!
                    for (int k=0;k<balls[i]._rimpos.size();k++)
                    {
                        if (balls[i]._times[k]>t)
                        {
                            break;
                        }
                        //check for each point here.
                        if (j==1)
                        {
                            a2=sqrt(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                            a3=sqrt(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                        }
                        else
                        {
                            a2=sqrt(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                            a3=sqrt(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                        }
                        if (((a2<ball_radius+2.094 && b2>0.) || (a3<ball_radius+2.094 && b3>0.)))
                        {
                            t=balls[i]._times[k];
                            collision=true;
                            //std::cout << "Hit" << std::endl;
                            break;
                        }
                    }
                    c=j+1;
                    break;
                }
            }
            //std::cout << "Checked middle pockets" << std::endl;
            if (c==0)
            {
                for (int j=0;j<4;j++)
                {
                    a1=sqrt(pow(balls[i]._x-cpockets[j][0],2.)+pow(balls[i]._y-cpockets[j][1],2.));
                    if (a1<cpocket_r+epsilon)
                    {
                        //rolling in this pocket!
                        for (int k=0;k<balls[i]._rimpos.size();k++)
                        {
                            if (balls[i]._times[k]>t)
                            {
                                break;
                            }
                            //check for each point here.
                            //check curved bit first.
                            if (j==0)
                            {
                                a2=sqrt(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==1)
                            {
                                a2=sqrt(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==2)
                            {
                                a2=sqrt(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==3)
                            {
                                a2=sqrt(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            if (((a2<ball_radius+4.5 && b2>0.) || (a3<ball_radius+4.5 && b3>0.)))
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                            //check straight bits.
                            t0=0.5*cpline[2*j][0]*(balls[i]._rimpos[k][1]+cpline[2*j][0]*balls[i]._rimpos[k][0]-cpline[2*j][1]);
                            t1=cpline[2*j][0]*t0+cpline[2*j][1];
                            a2=sqrt(pow(balls[i]._rimpos[k][0]-t0,2.)+pow(balls[i]._rimpos[k][1]-t1,2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t0-balls[i]._rimpos[k][0],t1-balls[i]._rimpos[k][1]));

                            t2=0.5*cpline[2*j+1][0]*(balls[i]._rimpos[k][1]+cpline[2*j+1][0]*balls[i]._rimpos[k][0]-cpline[2*j+1][1]);
                            t3=cpline[2*j+1][0]*t2+cpline[2*j+1][1];
                            a3=sqrt(pow(balls[i]._rimpos[k][0]-t2,2.)+pow(balls[i]._rimpos[k][1]-t3,2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t2-balls[i]._rimpos[k][0],t3-balls[i]._rimpos[k][1]));
                            if (((a2<ball_radius && t0>cpline[2*j][2] && t0<cpline[2*j][3] && b2>0.) || (a3<ball_radius && t2>cpline[2*j+1][2] && t2<cpline[2*j+1][3] && b3>0.)))
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                        }
                        c=j+3;
                        break;
                    }
                }
            }

            //check for collision with other balls.
            for (int j=0;j<22;j++)
            {
                if (i==j || balls[j]._potted)
                {
                    continue;
                }
                if (balls[j]._onrim)
                {
                    //check which pocket.
                    for (int k=0;k<2;k++)
                    {
                        a1=sqrt(pow(balls[j]._x-mpockets[k][0],2.)+pow(balls[j]._y-mpockets[k][1],2.));
                        if (a1<mpocket_r)
                        {
                            if (c==k+1)
                            {
                                //same pocket!
                                //check if and when they collide.
                                for (int m=0;m<std::min(balls[i]._times.size(),balls[j]._times.size());m++)
                                {
                                    if (balls[i]._times[m]>t)
                                    {
                                        break;
                                    }
                                    a1=sqrt(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2.)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2.));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2.*ball_radius && a2>0.)
                                    {
                                        t=balls[i]._times[m];
                                        collision=true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    for (int k=0;k<4;k++)
                    {
                        a1=sqrt(pow(balls[j]._x-cpockets[k][0],2.)+pow(balls[j]._y-cpockets[k][1],2.));
                        if (a1<cpocket_r)
                        {
                            if (c==k+3)
                            {
                                //same pocket!
                                //check if and when they collide.
                                for (int m=0;m<std::min(balls[i]._times.size(),balls[j]._times.size());m++)
                                {
                                    if (balls[i]._times[m]>t)
                                    {
                                        break;
                                    }
                                    a1=sqrt(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2.)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2.));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2.*ball_radius && a2>0.)
                                    {
                                        t=balls[i]._times[m];
                                        collision=true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                else
                {
                    xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                    xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                    xmin2-=1.01*ball_radius;
                    xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                    xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                    xmax2+=1.01*ball_radius;
                    ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                    ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                    ymin2-=1.01*ball_radius;
                    ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                    ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                    ymax2+=1.01*ball_radius;
                    //check bounds.
                    if (c<3)
                    {
                        //middle pocket.
                        if (xmin2>blue_x+1.96 || xmax2<blue_x-1.96)
                        {
                            continue;
                        }
                        if ((c==1 && ymax2<mpockets[0][1]-mpocket_r) || (c==2 && ymin2>mpockets[1][1]+mpocket_r))
                        {
                            continue;
                        }
                    }
                    else
                    {
                        //corner pocket.
                        if (c<5 && xmin2>rail_thickness+5.)
                        {
                            continue;
                        }
                        else if (c>4 && xmax2<rail_thickness+2.*cush_thickness+table_length-5.)
                        {
                            continue;
                        }
                        if ((c==3 || c==5) && ymin2>rail_thickness+5.)
                        {
                            continue;
                        }
                        else if ((c==4 || c==6) && ymax2<rail_thickness+2.*cush_thickness+table_width-5.)
                        {
                            continue;
                        }
                    }

                    //potential collision.
                    for (int k=0;k<balls[i]._rimpos.size();k++)
                    {
                        if (balls[i]._times[k]>t)
                        {
                            break;
                        }
                        //check for distance.
                        //assume a 2d collision to approximate.
                        out[0]=balls[j]._ax[2]*balls[i]._times[k]*balls[i]._times[k]+balls[j]._ax[1]*balls[i]._times[k]+balls[j]._ax[0];
                        out[1]=balls[j]._ay[2]*balls[i]._times[k]*balls[i]._times[k]+balls[j]._ay[1]*balls[i]._times[k]+balls[j]._ay[0];
                        if (sqrt(pow(balls[i]._rimpos[k][0]-out[0],2.)+pow(balls[i]._rimpos[k][1]-out[1],2.))-2.*ball_radius<epsilon)
                        {
                            //distance correct.
                            out2[0]=balls[j]._avx[1]*balls[i]._times[k]+balls[j]._avx[0];
                            out2[1]=balls[j]._avy[1]*balls[i]._times[k]+balls[j]._avy[0];
                            b1=atan2(balls[i]._rimpos[k][0]-out[0],balls[i]._rimpos[k][1]-out[1]);
                            a1=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b1)-sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-b1);
                            if (a1>0.)
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                        }
                    }
                }
            }
        }

        //add positions to list for each timestep.
        start=ceil(totaltime/0.01);
        for (int i=0;i<int(floor((totaltime+t)/0.01)-start)+1;i++)
        {
            t0=(start+i)*0.01-totaltime;
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted)
                {
                    temp[3*j]=-100.;
                    temp[3*j+1]=-100.;
                    temp[3*j+2]=-100.;
                    continue;
                }
                if (!balls[j]._onrim)
                {
                    temp[3*j]=balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0];
                    temp[3*j+1]=balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0];
                    temp[3*j+2]=balls[j]._az[2]*t0*t0+balls[j]._az[1]*t0+balls[j]._az[0];
                }
                else
                {
                    if (int(floor(1000.*t0))==balls[i]._rimpos.size()-1)
                    {
                        out=balls[i]._rimpos[balls[i]._rimpos.size()-1];
                        temp[3*j]=out[0];
                        temp[3*j+1]=out[1];
                        temp[3*j+2]=out[2];
                    }
                    else
                    {
                        //interpolate.
                        a1=balls[j]._times[int(floor(1000.*t0)+1)]-balls[j]._times[int(floor(1000.*t0))];
                        out=balls[j]._rimpos[int(floor(1000.*t0))];
                        out2=balls[j]._rimpos[int(floor(1000.*t0)+1)];
                        temp[3*j]=out[0]+(out2[0]-out[0])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
                        temp[3*j+1]=out[1]+(out2[1]-out[1])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
                        temp[3*j+2]=out[2]+(out2[2]-out[2])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
                    }
                }
            }
            pos.push_back(temp);
        }
        totaltime+=t;

        //fast forward to the next event. new velocity and displacement for the balls.
        //make sure the balls are positioned to collide properly (rounding error?)
        //update spins here.

        grid={};
        grid_index={};
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted)
            {
                continue;
            }
            if (!balls[i]._onrim)
            {
                //position.
                balls[i]._x=balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0];
                balls[i]._y=balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0];
                balls[i]._z=balls[i]._az[2]*t*t+balls[i]._az[1]*t+balls[i]._az[0];
                if (balls[i]._z<-ball_radius)
                {
                    balls[i]._potted=true;
                    balls[i]._vx=0.;
                    balls[i]._vy=0.;
                    balls[i]._vz=0.;
                    balls[i]._xspin=0.;
                    balls[i]._yspin=0.;
                    balls[i]._rspin=0.;
                    pos[pos.size()-1][3*i]=-100.;
                    pos[pos.size()-1][3*i+1]=-100.;
                    pos[pos.size()-1][3*i+2]=-100.;
                    ball_potted_order.push_back(balls[i]._order);
                }
                if (!balls[i]._potted)
                {
                    //velocity.
                    balls[i]._vx=balls[i]._avx[1]*t+balls[i]._avx[0];
                    balls[i]._vy=balls[i]._avy[1]*t+balls[i]._avy[0];
                    balls[i]._vz=balls[i]._avz[1]*t+balls[i]._avz[0];
                    if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                    if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                    if (fabs(balls[i]._vz)<epsilon) {balls[i]._vz=0.;}
                    //spin.
                    balls[i]._xspin=balls[i]._awx[1]*t+balls[i]._awx[0];
                    balls[i]._yspin=balls[i]._awy[1]*t+balls[i]._awy[0];
                    balls[i]._rspin=balls[i]._awz[1]*t+balls[i]._awz[0];
                    if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0.;}
                    if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0.;}
                    if (fabs(balls[i]._rspin)<epsilon) {balls[i]._rspin=0.;}
                }
            }
            else
            {
                //on the rim.
                if (int(floor(1000.*t))==balls[i]._rimpos.size()-1)
                {
                    out=balls[i]._rimpos[balls[i]._rimpos.size()-1];
                    balls[i]._x=out[0];
                    balls[i]._y=out[1];
                    balls[i]._z=out[2];
                    out=balls[i]._rimvel[balls[i]._rimvel.size()-1];
                    balls[i]._vx=out[0];
                    balls[i]._vy=out[1];
                    balls[i]._vz=out[2];
                }
                else
                {
                    //interpolate.
                    a1=balls[i]._times[int(floor(1000.*t)+1)]-balls[i]._times[int(floor(1000.*t))];
                    out=balls[i]._rimpos[int(floor(1000.*t))];
                    out2=balls[i]._rimpos[int(floor(1000.*t)+1)];
                    balls[i]._x=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._y=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._z=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    out=balls[i]._rimvel[int(floor(1000.*t))];
                    out2=balls[i]._rimvel[int(floor(1000.*t)+1)];
                    balls[i]._vx=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._vy=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._vz=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                }
                balls[i]._xspin=balls[i]._vx/ball_radius;
                balls[i]._yspin=balls[i]._vy/ball_radius;
                //work out how rspin changes.
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                if (fabs(balls[i]._vz)<epsilon) {balls[i]._vz=0.;}
                if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0.;}
                if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0.;}
            }
            //update the grid positions in case there is a collision.
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
            //update grid position.

            for (int j=0;j<4;j++)
            {
                xpos=balls[i]._gpos[j][0];
                ypos=balls[i]._gpos[j][1];
                grid[xpos][ypos][grid_index[xpos][ypos]]=balls[i]._order;
                grid_index[xpos][ypos]+=1;
            }

//            if (i==5)
//            {
//                std::cout << "x:" << balls[i]._x << std::endl;
//                std::cout << "y:" << balls[i]._y << std::endl;
//                std::cout << "z:" << balls[i]._z << std::endl;
//                std::cout << "vx:" << balls[i]._vx << std::endl;
//                std::cout << "vy:" << balls[i]._vy << std::endl;
//                std::cout << "vz:" << balls[i]._vz << std::endl;
//                std::cout << "xspin:" << double(balls[i]._xspin) << std::endl;
//                std::cout << "yspin:" << double(balls[i]._yspin) << std::endl;
//                std::cout << "rspin:" << double(balls[i]._rspin) << std::endl;
//            }
        }

        //perform collision check.
        if (collision)
        {
            do
            {
                dv=collisions(balls,cush);
                for (int i=0;i<22;i++)
                {
                    if (balls[i]._potted)
                    {
                        continue;
                    }
                    balls[i]._vx+=dv(2*i);
                    balls[i]._vy+=dv(2*i+1);
                    if (fabs(dv(2*i))>=epsilon && fabs(dv(2*i+1))>=epsilon && balls[i]._inflight)
                    {
                        a1=999.;
                        a2=999.;
                        for (int k=0;k<4;k++)
                        {
                            a1=fmin(a1,sqrt(pow(balls[i]._x-cpockets[k][0],2.)+pow(balls[i]._y-cpockets[k][1],2.)));
                        }
                        for (int k=0;k<2;k++)
                        {
                            a2=fmin(a2,sqrt(pow(balls[i]._x-mpockets[k][0],2.)+pow(balls[i]._y-mpockets[k][1],2.)));
                        }
                        if (a1<k1 || (a2<mpocket_r && fabs(blue_y-balls[i]._y)>0.5*table_width+0.438-ball_radius-0.01))
                        {
                            balls[i]._vz=-2*sqrt(pow(balls[i]._vx,2.)+pow(balls[i]._vy,2.));
                        }
                    }
                }
                if (!hitdone)
                {
                    if (fabs(dv(0))>=epsilon || fabs(dv(1))>=epsilon)
                    {
                        for (int j=1;j<22;j++)
                        {
                            if (fabs(dv(2*j))>=epsilon || fabs(dv(2*j+1))>=epsilon)
                            {
                                ball_hit_order.push_back(balls[j]._order);
                            }
                        }
                        hitdone=true;
                    }
                }
            } while (dv!=zero);
            for (int i=0;i<22;i++)
            {
                if (balls[i]._potted)
                {
                    continue;
                }
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                for (int j=i+1;j<22;j++)
                {
                    t0=sqrt(pow(balls[i]._x-balls[j]._x,2.)+pow(balls[i]._y-balls[j]._y,2.));
                    if (t0-2.*ball_radius<=0.)
                    {
                        t1=atan2(balls[j]._x-balls[i]._x,balls[j]._y-balls[i]._y);
                        t2=(2.*ball_radius-t0)+DOUBLE_EPSILON;
                        balls[i]._x-=t2*sin(t1);
                        balls[i]._y-=t2*cos(t1);
                        balls[j]._x+=t2*sin(t1);
                        balls[j]._y+=t2*cos(t1);
                    }
                }
            }
        }
    }
    return pos;
}

void Server::handleIncomingConnections()
{
    bool alreadyconnected=false;
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()==sf::IpAddress::None)
        {
            //fresh socket ready for connecting to.
            if (listener.accept(players[i])==sf::TcpListener::Done)
            {
                //connected to new socket.
                alreadyconnected=true;
                break;
            }
        }
    }
    if (!alreadyconnected)
    {
        for (int i=0;i<4;i++)
        {
            if (spectators[i].getRemoteAddress()==sf::IpAddress::None)
            {
                //fresh socket ready for connecting to.
                if (listener.accept(spectators[i])==sf::TcpListener::Done)
                {
                    //connected to new socket.
                    break;
                }
            }
        }
    }
}

void Server::executionThread()
{
    double a=0.;
    double b=0.;
    double c=0.;
    double d=0.;
    double e=0.;
    double f=0.;
    double g=0.;
    std::string msg;

    while (running)
    {
        handleIncomingConnections();

        //check for incoming packets.
        //only receive input from the active player.
        if (players[player_turn].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            if (players[player_turn].receive(packet)==sf::Socket::Done)
            {
                //received an input packet.
                packet >> packetId;
                if (packetId==0)
                {
                    //disconnect this user.
                }
                else if (packetId==1)
                {
                    //trajectory display to other clients.
                    packet >> a >> b >> c >> d >> e >> f >> g >> nom_colour_order;
                    packet.clear();
                    packet << sf::Uint16(0) << a << b << c << d << e << f << g << nom_colour_order;

                    serverballs[0]._x=f;
                    serverballs[0]._y=g;

                    if (players[!player_turn].getRemoteAddress()!=sf::IpAddress::None)
                    {
                        players[!player_turn].send(packet);
                    }

                    for (int i=0;i<4;i++)
                    {
                        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
                        {
                            spectators[i].send(packet);
                        }
                    }
                }
                else if (packetId==2)
                {
                    //store original position.
                    for (int i=0;i<22;i++)
                    {
                        potbefore[i]=serverballs[i]._potted;
                        posbefore[3*i]=serverballs[i]._x;
                        posbefore[3*i+1]=serverballs[i]._y;
                        posbefore[3*i+2]=serverballs[i]._z;
                    }
                    wasredon=isredon;
                    wasfreeball=isfreeball;
                    wasRedsLeft=redsLeft;
                    //simulate and return the results to everybody.
                    packet >> a >> b >> c >> d >> e;
                    serverballs[0]._vx=a;
                    serverballs[0]._vy=b;
                    serverballs[0]._xspin=c;
                    serverballs[0]._yspin=d;
                    serverballs[0]._rspin=e;
                    result.clear();
                    result=simulate(serverballs,servercushions);

                    respot();

                    applyRulesAndScore();

                    //send the log message if foul or miss.
                    if (isfoul)
                    {
                        if (ismiss)
                        {
                            msg="FOUL AND A MISS - "+pnames[!player_turn]+", enter /1 to play from here, /2 for opponent to play from here, or /3 for opponent to play from original position.";
                            broadcast("~",msg);
                        }
                        else
                        {
                            msg="FOUL - "+pnames[!player_turn]+", enter /1 to play from here or /2 for opponent to play from here.";
                            broadcast("~",msg);
                        }

                        if (isfreeball)
                        {
                            broadcast("~","Free ball (if play from here).");
                        }
                    }

                    sendShotSimulation();

                    turnpacket();

                    //check if all balls have been potted, then reset the frame.
                    if (!redsLeft)
                    {
                        bool allPotted=true;
                        for (int i=1;i<22;i++)
                        {
                            if (serverballs[i]._potted==false) {allPotted=false; break;}
                        }
                        if (allPotted==true)
                        {
                            resetframe();
                            turnpacket();
                        }
                    }
                }
                else if (packetId==3)
                {
                    //message to the server.
                    packet >> msg;

                    broadcast(pnames[player_turn],msg);

                    if (isfoul)
                    {
                        handleFoulMissResponse(msg);
                    }
                }
                else if (packetId==4)
                {
                    //concede frame.
                    frames[!player_turn]+=1;
                    if (frames[!player_turn]==(framesbestof+1)/2)
                    {
                        gameover=true;
                    }
                    resetframe();
                    turnpacket();
                }
                else if (packetId==5)
                {
                    //concede match.
                    frames[!player_turn]=(framesbestof+1)/2;
                    gameover=true;
                    turnpacket();
                }
                else if (packetId==6)
                {
                    //name of player.
                    packet >> msg;
                    pnames[player_turn]=msg;
                    broadcast("~","Welcome, "+pnames[player_turn]+"!");
                    sendNamePacket();
                    sendBallPositions();
                    turnpacket();
                }
                else if (packetId==7)
                {
                    //update placing_white.
                    packet >> placing_white;
                    turnpacket();
                }
            }
        }

        if (players[!player_turn].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            if (players[!player_turn].receive(packet)==sf::Socket::Done)
            {
                packet >> packetId;
                if (packetId==3)
                {
                    packet >> msg;

                    broadcast(pnames[!player_turn],msg);
                }
                else if (packetId==4)
                {
                    //concede frame.
                    frames[player_turn]+=1;
                    if (frames[player_turn]==(framesbestof+1)/2)
                    {
                        gameover=true;
                    }
                    resetframe();
                    turnpacket();
                }
                else if (packetId==5)
                {
                    //concede match.
                    frames[player_turn]=(framesbestof+1)/2;
                    gameover=true;
                    turnpacket();
                }
                else if (packetId==6)
                {
                    //name of player.
                    packet >> msg;
                    pnames[!player_turn]=msg;
                    broadcast("~","Welcome, "+pnames[!player_turn]+"!");
                    sendNamePacket();
                    sendBallPositions();
                    turnpacket();
                }
            }
        }

//        for (int i=0;i<4;i++)
//        {
//            if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
//            {
//                packet.clear();
//                if (spectators[i].receive(packet)==sf::Socket::Done)
//                {
//                    packet >> packetId;
//                    if (packetId==3)
//                    {
//                        packet >> msg;
//
//                        broadcast(snames[i],msg);
//                    }
//                }
//            }
//        }

        sf::sleep(sf::milliseconds(5));
    }
}

void Server::resetServer()
{
    for (int i=0;i<players.size();i++)
    {
        players[i].disconnect();
    }
    for (int i=0;i<spectators.size();i++)
    {
        spectators[i].disconnect();
    }
    player_turn=0;
    pnames[0]="PLAYER 1";
    pnames[1]="PLAYER 2";
    scores[0]=0;
    scores[1]=0;
    frames[0]=0;
    frames[1]=0;
    framesbestof=3;
    placing_white=true;
    touching=false;
    isfoul=false;
    isredon=true;
    nom_colour_order=0; //1 yellow,2 green, etc.
    isfreeball=false;

    gameover=false;
    highbreak[0]=0;
    highbreak[1]=0;
    centuries[0]=0;
    centuries[1]=0;

    rackballs();
}

#endif // SERVER_H_INCLUDED
