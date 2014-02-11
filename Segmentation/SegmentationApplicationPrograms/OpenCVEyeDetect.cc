//File: OpenCVFaceDetect.cc
//Input: Query Image Filename, Cascade Path
//Output: Displays the input image with a rectangle around it to show the detected face
//Author: Bud Goswami
//Date: 26.01.09
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/String.hh"
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include "Ravl/Point2d.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Image/DrawCross.hh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>

using namespace RavlImageN;

//OPENCV STUFF
// Create memory for calculations
//static CvMemStorage* storage = 0;
// Create a new Haar classifier
CvHaarClassifierCascade* cascade = 0;
CvHaarClassifierCascade* cascade_e = 0;
// Function prototype for detecting and drawing an object from an image
ImageRectangleC detect_and_draw( IplImage* img );
Tuple2C<Point2dC,Point2dC> detectEyes(IplImage *img);
//NORMAL RAVL STUFF

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory");
	FilenameC cascade_name = opt.String("c","/opt/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	FilenameC eye_cascade_name = opt.String("e","/vol/vssp/lip-tracking/OpenCVFaceDetect/Eyes/haarcascade_eye.xml","Path for the Haar Classifier");
	IntT display = opt.Int("d",0, "Set to 1 if you want to save the images as well");
	opt.Check();
	
	DListC<StringC> qfile = qimg_dir.FiltList("*.ppm");
	for(DLIterC<StringC> it(qfile); it; it++)
	{
		FilenameC qimg = qimg_dir+(*it);
		
		ImageC<RealRGBValueC> src;   
		if(!Load(qimg, src)) cerr<<"Loading RAVL Image Failed"<<endl;
		IplImage* image = 0;
		image = cvLoadImage( qimg.chars(), 1 );
		cascade = (CvHaarClassifierCascade*)cvLoad( cascade_name.chars(), 0, 0, 0 ); 
		cascade_e = (CvHaarClassifierCascade*)cvLoad( eye_cascade_name.chars(), 0, 0, 0 ); 
		// Check whether the cascade has loaded successfully. Else report and error and quit
		if( !cascade )
		{
			fprintf( stderr, "ERROR: Could not load classifier cascade\n" );
			return -1;
		} 
		//CvMemStorage* storage = cvCreateMemStorage(0);
		ImageRectangleC facerec;Tuple2C<Point2dC, Point2dC> eyes;
		if( image )
		{
			//Detect and draw the face
			facerec = detect_and_draw( image );
			eyes = detectEyes( image );
			// Release the image memory
			cvReleaseImage( &image );
		}
		//eyes.Data1() = Point2dC(eyes.Data1().Row() + (RealT)facerec.TopLeft().Row(),eyes.Data1().Col() + (RealT)facerec.TopLeft().Col());
		cout<<"Point = "<< Point2dC(eyes.Data1().Row() + (RealT)facerec.TopLeft().Row(),eyes.Data1().Col() + (RealT)facerec.TopLeft().Col()) << endl;
		bool Fill = false;
		RealRGBValueC border(100,10,10);
		DrawFrame(src,border,facerec,Fill);
		DrawCross(src,border,eyes.Data1() + facerec.TopLeft(),3);    
		DrawCross(src,border,eyes.Data2() + facerec.TopLeft(),5);
		if(!Save("@X: Face with image rectangle",src)) cerr<<"Could not show face border"<<endl;
		cout<<"Enter number to continue"<<endl;IntT d;
		cin>>d;
		//	if(!Save(outimg, out)) cerr<<"Failed to save output image"<<endl;
		//	if(!Save(gtoutimg, outgt)) cerr<<"Failed to save output gt image"<<endl;		
	}
	return 0;
}

// Function prototype for detecting and drawing an object from an image
ImageRectangleC detect_and_draw( IplImage* img )
{
	Point2dC topleft(0,0); Point2dC botright(0,0);
	int scale = 1;
	CvMemStorage* storage = cvCreateMemStorage(0);
	 // Create a new image based on the input image
	 IplImage* temp = cvCreateImage( cvSize(img->width/scale,img->height/scale), 8, 3 );

	 // Create two points to represent the face locations
	 CvPoint pt1, pt2;
	 int i;
	 // Clear the memory storage which was used before
	 cvClearMemStorage(storage);
	 // Find whether the cascade is loaded, to find the faces. If yes, then:
	 if( cascade )
	 {
		 	// There can be more than one face in an image. So create a growable sequence of faces.
	        // Detect the objects and store them in the sequence
	        CvSeq* faces = cvHaarDetectObjects( img, cascade, storage,1.1, 2, CV_HAAR_DO_CANNY_PRUNING,cvSize(40, 40) );
	        // Loop the number of faces found.
	        for( i = 0; i < (faces ? faces->total : 0); i++ )
	        {
	        	// Create a new rectangle for drawing the face
	            CvRect* r = (CvRect*)cvGetSeqElem( faces, i );
	            // Find the dimensions of the face,and scale it if necessary
	            pt1.x = r->x*scale;
	            pt2.x = (r->x+r->width)*scale;
	            pt1.y = r->y*scale;
	            pt2.y = (r->y+r->height)*scale;
	            topleft.Row() = pt1.y;
	            topleft.Col() = pt1.x;
	            botright.Row() = pt2.y;
	            botright.Col() = pt2.x;            
	            // Draw the rectangle in the input image
	            cvRectangle( img, pt1, pt2, CV_RGB(255,0,0), 3, 8, 0 );
	        }
	    }
	    // Release the temp image created.
	    cvReleaseImage( &temp );
		return ImageRectangleC(topleft,botright);
}

Tuple2C<Point2dC, Point2dC> detectEyes(IplImage *img)
{
	Point2dC l_eye(0.0,0.0),r_eye(0.0,0.0);
	CvMemStorage* storage = cvCreateMemStorage(0);
	int i;
    /* detect faces */
	CvSeq *faces = cvHaarDetectObjects(img, cascade, storage,1.1, 3, 0, cvSize( 40, 40 ) );
    /* return if not found */
    //if (faces->total == 0) return;
    /* draw a rectangle */
	CvRect *r = (CvRect*)cvGetSeqElem(faces, 0);
	cvRectangle(img, cvPoint(r->x, r->y),cvPoint(r->x + r->width, r->y + r->height),CV_RGB(255, 0, 0), 1, 8, 0);
    /* reset buffer for the next object detection */
    cvClearMemStorage(storage);
    /* Set the Region of Interest: estimate the eyes' position */
    cvSetImageROI(img, cvRect(r->x, r->y + (r->height/5.5), r->width, r->height/3.0));
    /* detect eyes */
	CvSeq* eyes = cvHaarDetectObjects(img, cascade_e, storage,1.15, 3, 0, cvSize(25, 15));
	r = (CvRect*)cvGetSeqElem( eyes, 0 ); //get left eye
	l_eye.Row() = r->y; l_eye.Col() = r->x;
	r = (CvRect*)cvGetSeqElem( eyes, 0 ); //get right eye
	r_eye.Row() = r->y; r_eye.Col() = r->x;
	Tuple2C<Point2dC, Point2dC> res(l_eye.Copy(),r_eye.Copy());
	cout<<"Found eyes = "<<res<<endl;
	return res;
	/* draw a rectangle for each eye found */
	/*
	for( i = 0; i < (eyes ? eyes->total : 0); i++ ) {
		r = (CvRect*)cvGetSeqElem( eyes, i );
		cvRectangle(img, cvPoint(r->x, r->y), cvPoint(r->x + r->width, r->y + r->height),CV_RGB(255, 0, 0), 1, 8, 0);
	}
    cvResetImageROI(img);
	*/
}
