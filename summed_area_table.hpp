#ifndef SUMTABLE_H_
#define SUMTABLE_H_

#include <imutil.hpp>

using namespace imutil;

template<typename TYPE_IN, typename TYPE_OUT>
class summed_area_table
{
	private:
        void build(const cv::Mat& in) {
            table=cv::Mat::zeros(in.size(),cv::DataType<TYPE_OUT>::type);

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

        summed_area_table(const cv::Mat& in) {
            build(in);
        }

        TYPE_OUT sum(int x,int y,int r) {
            int xbr=x+r,ybr=y+r;
            int xtl=x-r-1,ytl=y-r-1;
            int xtr=x+r,ytr=y-r-1;
            int xbl=x-r-1,ybl=y+r;

            if(xbr>=table.cols) xbr=table.cols-1;
            if(ybr>=table.rows) ybr=table.rows-1;

            TYPE_OUT br=imutil::get_px<TYPE_OUT>(table,xbr,ybr);
            TYPE_OUT winsum=br;

            if(xtl>=0&&ytl>=0) {
                TYPE_OUT tl=imutil::get_px<TYPE_OUT>(table,xtl,ytl);
                winsum+=tl;
            }

            if(xbl>=0) {
                if(ybl>=table.rows) ybl=table.rows-1;

                TYPE_OUT bl=imutil::get_px<TYPE_OUT>(table,xbl,ybl);
                winsum-=bl;
            }

            if(ytr>=0) {
                if(xtr>=table.cols) xtr=table.cols-1;

                TYPE_OUT tr=imutil::get_px<TYPE_OUT>(table,xtr,ytr);
                winsum-=tr;
            }

            return winsum;
        }

        TYPE_OUT mean(int x,int y,int r) {
            int xbr=x+r,ybr=y+r;
            int xtl=x-r-1,ytl=y-r-1;

            if(xbr>=table.cols) xbr=table.cols-1;
            if(ybr>=table.rows) ybr=table.rows-1;

            int wtlx=xtl>=-1?xtl:-1,wtly=ytl>=-1?ytl:-1,wbrx=xbr,wbry=ybr;
            int npx=(wbrx-wtlx)*(wbry-wtly);
            TYPE_OUT winsum=sum(x,y,r);

            return winsum/npx;
        };

        //when window size is large (>80) you get strange results (the variance is always 0 etc.)
        TYPE_OUT variance(summed_area_table& table_squared,int x,int y,int r) {
            TYPE_OUT S1=sum(x,y,r);
            TYPE_OUT S2=table_squared.sum(x,y,r);
            int xbr=x+r,ybr=y+r;
            int xtl=x-r-1,ytl=y-r-1;

            if(xbr>=table.cols) xbr=table.cols-1;
            if(ybr>=table.rows) ybr=table.rows-1;

            int wtlx=xtl>=-1?xtl:-1,wtly=ytl>=-1?ytl:-1,wbrx=xbr,wbry=ybr;
            int npx=(wbrx-wtlx)*(wbry-wtly);

            return (S2-pow(S1,2)/npx)/npx;
        };
};

#endif
