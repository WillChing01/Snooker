#ifndef CUSTOM-RECEIVER_H_INCLUDED
#define CUSTOM-RECEIVER_H_INCLUDED

class CustomReceiver: public sf::SoundStream
{
    private:
        unsigned short port;
        sf::IpAddress ipAddress;
        sf::TcpSocket socket; //blocking by default.
        sf::Packet packet;

        virtual bool onGetData(Chunk& data)
        {
            //data.samples =
            //data.sampleCount =

            return true; //return false to stop playback.
        }
        virtual void onSeek(sf::Time timeOffset)
        {
            //do nothing yet.
        }
    public:
        void setTarget(sf::IpAddress newIp, unsigned short newPort)
        {
            ipAddress=newIp;
            port=newPort;
            try
            {
                socket.connect(ipAddress,port);
            }
            catch (...) {}
        }
};

#endif // CUSTOM-RECEIVER_H_INCLUDED
