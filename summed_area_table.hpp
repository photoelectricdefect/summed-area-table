#ifndef SUMTABLE_H_
#define SUMTABLE_H_

#include <imutil.hpp>

template<typename TYPE_IN, typename TYPE_OUT>
class summed_area_table
{
	private:
		void build(const cv::Mat& in) {
			table=cv::Mat::zeros(in.size(),TYPE_OUT);

			for(size_t x=0;x<in.cols;x++) {
            	TYPE_OUT sum=get_px<TYPE_IN>(in,x,0);

            	if(x>0) sum+=get_px<TYPE_OUT>(table,x-1,0);

            	set_px<TYPE_OUT>(table,x,0,sum);
        	}

        	for(size_t y=0;y<in.rows;y++) {
            	TYPE_OUT sum=get_px<TYPE_IN>(in,0,y);

            	if(y>0) sum+=get_px<TYPE_OUT>(table,0,y-1);  

            	set_px<TYPE_OUT>(table,0,y,sum);
        	}

        	for(size_t x=1;x<in.cols;x++) {
            	for(size_t y=1;y<in.rows;y++) {
                	TYPE_OUT C =get_px<TYPE_IN>(in,x,y),Ix=get_px<TYPE_OUT>(table,x-1,y),Iy=get_px<TYPE_OUT>(table,x,y-1),Ixy=get_px<TYPE_OUT>(table,x-1,y-1);
                	TYPE_OUT I=C+Ix+Iy-Ixy;
                	set_px<TYPE_OUT>(table,x,y,I);            
            	}
        	}
    	}
	public:
		cv::Mat table;
		summed_area_table(const cv::Mat& in) {build(in);}
};

#endif