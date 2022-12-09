#ifndef AUDIO-SERVER_H_INCLUDED
#define AUDIO-SERVER_H_INCLUDED

class AudioServer
{
    public:
        bool running=true;

        sf::IpAddress serverIp;
        unsigned short port;
        sf::TcpListener listener;

        sf::Packet packet;
        sf::Uint16 packetId=0;

        std::array<sf::TcpSocket,6> clients;

        AudioServer(bool isOnline=true);
        void broadcast(sf::Packet packet);
        void handleIncomingConnections();
        void executionThread();
};

AudioServer::AudioServer(bool isOnline)
{
    if (isOnline)
    {
        try {serverIp=sf::IpAddress::getLocalAddress();}
        catch (...) {}
        listener.setBlocking(false);
        if (listener.listen(sf::Socket::AnyPort)!=sf::Socket::Done)
        {
            //error.
            std::cout << "Tcp listener could not bind to port." << std::endl;
        }
        try{port=listener.getLocalPort();}
        catch (...) {}
    }

    for (int i=0;i<clients.size();i++)
    {
        clients[i].setBlocking(false);
    }
}

void AudioServer::broadcast(sf::Packet packet)
{
    for (int i=0;i<clients.size();i++)
    {
        clients[i].send(packet);
    }
}

void AudioServer::handleIncomingConnections()
{
    for (int i=0;i<clients.size();i++)
    {
        if (clients[i].getRemoteAddress()==sf::IpAddress::None)
        {
            //fresh socket ready for connecting to.
            if (listener.accept(clients[i])==sf::TcpListener::Done)
            {
                //connected to new socket.
                break;
            }
        }
    }
}

void AudioServer::executionThread()
{
    while (running)
    {
        handleIncomingConnections();

        for (int i=0;i<clients.size();i++)
        {
            if (clients[i].getRemoteAddress()==sf::IpAddress::None) {continue;}

            //check for incoming packets.
            if (clients[i].receive(packet)==sf::Socket::Done)
            {
                //received an audio packet.
                broadcast(packet);
            }
        }
    }
}

#endif // AUDIO-SERVER_H_INCLUDED
