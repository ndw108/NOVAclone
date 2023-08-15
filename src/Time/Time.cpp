/*  This file is part of NOVA.

    NOVA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NOVA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NOVA.  If not, see <http://www.gnu.org/licenses/>.

    Author: Alex Skillen. alex.skillen@manchester.ac.uk

*/
#include "Time/Time.h"
#include <iomanip>

Time::Time
(
    scalar dt,
    scalar endTime,
    int writeEvery
)
:
    steps_(),
    dt_(dt),
    curTime_(0.0),
    endTime_(endTime),
    writeEvery_(writeEvery),
    timeStep_(0),
    fileNum_(0),
    writeEvery_timebased_(false),
    write_(false),
    log_(true)
{
    if( settings::restart() )
    {
        std::ifstream fin;
        std::stringstream fileName;
        fileName << "./data/restart.p" << std::setfill('0') << std::setw(4) << parallelCom::myProcNo() << ".case";
        fin.open(fileName.str());

        unsigned n;
        fin>>n;

        double t;
        for( unsigned i=0 ; i<n ; i++ )
        {
            fin >> t;
            steps_.push_back( t );
        }

        curTime_ = t;

        fin.close();

        fileNum_ = steps_.size();
    }
}

Time::Time
(
    scalar dt,
    scalar endTime,
    scalar writeEvery
)
: 
    Time(dt, endTime, 0)
{
    writeEvery_timebased_ = true;
    writeInterval_=writeEvery;
}


Time::Time
(
    std::string fileName
)
:
    Time(0,0,0)
{
    std::string line;
    std::ifstream fin;
    fin.open( fileName ); 

    while(std::getline(fin, line))
    {
        std::stringstream strProcessor(line);
        std::string buf;
    
        strProcessor>>buf;

        if( buf == std::string( "dt:" ) )
        {
            strProcessor>>dt_;
        }
        else if( buf == std::string( "writeT:" ) )
        {
            writeEvery_timebased_ = true;
            strProcessor>>writeInterval_;
        }
        else if( buf == std::string( "writeN:" ) )
        {
            writeEvery_timebased_ = false;
            strProcessor>>writeEvery_;
        }
        else if( buf == std::string( "end:" ) )
        {
            strProcessor>>endTime_;
        }
    }
}


void Time::operator++(int)
{
    if( writeEvery_timebased_ )
    {
        write_ = ( (curTime_+dt_) >= (fileNum_+1)*writeInterval_ );
    }
    else
    {
        write_ = ( timeStep_ % writeEvery_ == 0 );
    }    

    if( this->write() )
    {
        fileNum_++;
        steps_.push_back( curTime_ );
    }

    log_ = ( timeStep_ % 100 == 0 );

    curTime_+=dt_;
    timeStep_++;
} 
