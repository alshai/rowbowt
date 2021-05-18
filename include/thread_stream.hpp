#ifndef THREADSTREAM
#define THREADSTREAM

#include <iostream>
#include <sstream>
#include <mutex>

#define terr ThreadStream(std::cerr)
#define tout ThreadStream(std::cout)

/**
 * Thread-safe std::ostream class.
 *
 * Usage:
 *    tout << "Hello world!" << std::endl;
 *    terr << "Hello world!" << std::endl;
 */
class ThreadStream : public std::ostringstream
{
    public:
        ThreadStream(std::ostream& os) : os_(os)
        {
            // copyfmt causes odd problems with lost output
            // probably some specific flag
//            copyfmt(os);
            // copy whatever properties are relevant
            imbue(os.getloc());
            precision(os.precision());
            width(os.width());
            setf(std::ios::fixed, std::ios::floatfield);
        }

        ~ThreadStream()
        {
            std::lock_guard<std::mutex> guard(_mutex_threadstream);
            os_ << this->str();
        }

    private:
        static std::mutex _mutex_threadstream;
        std::ostream& os_;
};

std::mutex ThreadStream::_mutex_threadstream{};

#endif
