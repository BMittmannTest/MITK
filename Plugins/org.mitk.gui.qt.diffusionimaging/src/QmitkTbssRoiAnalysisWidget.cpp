/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center, 
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR 
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkTbssRoiAnalysisWidget.h"


#include <qlabel.h>
#include <qpen.h>
#include <qgroupbox.h>

#include <qwt_text.h>
#include <qwt_plot_grid.h>
#include <qwt_painter.h>
#include <qwt_legend.h>
#include <qwt_plot_marker.h>


#include <iostream>
#include <fstream>
#include <limits>


QmitkTbssRoiAnalysisWidget::QmitkTbssRoiAnalysisWidget( QWidget * parent )
  : QmitkPlotWidget(parent)
{
  m_PlotPicker = new QwtPlotPicker(m_Plot->canvas());
  m_PlotPicker->setSelectionFlags(QwtPicker::PointSelection | QwtPicker::ClickSelection | QwtPicker::DragSelection);
  m_PlotPicker->setTrackerMode(QwtPicker::ActiveOnly);
}

std::vector< std::vector<double> > QmitkTbssRoiAnalysisWidget::CalculateGroupProfiles(std::string preprocessed)
{
  MITK_INFO << "make profiles!";
  std::vector< std::vector<double> > profiles;

  //No results were preprocessed, so they must be calculated now.
  if(preprocessed == "")
  {
  // Iterate through the 4th dim (corresponding to subjects)
  // and create a profile for every subject
    int size = m_Projections->GetVectorLength();
    for(int s=0; s<size; s++)
    {      
      // Iterate trough the roi
      std::vector<double> profile;
      RoiType::iterator it;
      it = m_Roi.begin();
      while(it != m_Roi.end())
      {
        itk::Index<3> ix = *it;
        profile.push_back(m_Projections->GetPixel(ix).GetElement(s));
        it++;
      }
      int pSize = profile.size();
      profiles.push_back(profile);
    }
  }
  else{ // Use preprocessed results
    std::ifstream file(preprocessed.c_str());
    if(file.is_open())
    {
      std::string line;
      while(getline(file,line))
      {
        std::vector<std::string> tokens;
        Tokenize(line, tokens);
        std::vector<std::string>::iterator it;
        it = tokens.begin();
        std::vector< double > profile;
        while(it != tokens.end())
        {
          std::string s = *it;
          profile.push_back (atof( s.c_str() ) );
          ++it;
        }
        profiles.push_back(profile);
      }
    }
  }


  m_IndividualProfiles = profiles;
  // Calculate the averages

  // Here a check could be build in to check whether all profiles have
  // the same length, but this should normally be the case if the input
  // data were corrected with the TBSS Module.

  std::vector< std::vector<double> > groupProfiles;

  std::vector< std::pair<std::string, int> >::iterator it;
  it = m_Groups.begin();

  int c = 0; //the current profile number

  int nprof = profiles.size();

  while(it != m_Groups.end() && profiles.size() > 0)
  {
    std::pair<std::string, int> p = *it;
    int size = p.second;

    //initialize a vector of the right length with zeroes
    std::vector<double> averageProfile;
    for(int i=0; i<profiles.at(0).size(); i++)
    {
      averageProfile.push_back(0.0);
    }

    // Average the right number of profiles

    for(int i=0; i<size; i++)
    {
      for(int j=0; j<averageProfile.size(); ++j)
      {
        averageProfile.at(j) = averageProfile.at(j) + profiles.at(c).at(j);
      }
      c++;
    }

    // Devide by the number of profiles to get group average
    for(int i=0; i<averageProfile.size(); i++)
    {
      averageProfile.at(i) = averageProfile.at(i) / size;
    }

    groupProfiles.push_back(averageProfile);

    ++it;
  }
  int pSize = groupProfiles.size();

  return groupProfiles;
}


void QmitkTbssRoiAnalysisWidget::DrawProfiles(std::string preprocessed)
{
  this->Clear();



  m_Vals.clear();

  std::vector<double> v1;


  std::vector <std::vector<double> > groupProfiles = CalculateGroupProfiles(preprocessed);


  std::vector<double> xAxis;
  for(int i=0; i<groupProfiles.at(0).size(); ++i)
  {
    xAxis.push_back((double)i);
  }


  // fill m_Vals. This might be used by the user to copy data to the clipboard
  for(int i=0; i<groupProfiles.size(); i++)
  {

    v1.clear();


    for(int j=0; j<groupProfiles.at(i).size(); j++)
    {
      v1.push_back(groupProfiles.at(i).at(j));
    }

    m_Vals.push_back(v1);

  }

  std::string title = m_Measure + " profiles on the ";
  title.append(m_Structure);
  this->SetPlotTitle( title.c_str() );
  QPen pen( Qt::SolidLine );
  pen.setWidth(2);



  std::vector< std::pair<std::string, int> >::iterator it;
  it = m_Groups.begin();

  int c = 0; //the current profile number

  QColor colors[4] = {Qt::green, Qt::blue, Qt::yellow, Qt::red};

  while(it != m_Groups.end() && groupProfiles.size() > 0)
  {

    std::pair< std::string, int > group = *it;

    pen.setColor(colors[c]);
    int curveId = this->InsertCurve( group.first.c_str() );
    this->SetCurveData( curveId, xAxis, groupProfiles.at(c) );


    this->SetCurvePen( curveId, pen );

    c++;
    it++;

  }


  QwtLegend *legend = new QwtLegend;
  this->SetLegend(legend, QwtPlot::RightLegend, 0.5); 


  std::cout << m_Measure << std::endl;
  this->m_Plot->setAxisTitle(0, m_Measure.c_str());
  this->m_Plot->setAxisTitle(3, "Position");

  this->Replot();

}



void QmitkTbssRoiAnalysisWidget::PlotFiberBundles(TractContainerType tracts, mitk::Image *img)
{
  this->Clear();



  std::vector<TractType>::iterator it = tracts.begin();


  // Match points on tracts. Take the smallest tract and match all others on this one


  int min = std::numeric_limits<int>::max();
  TractType smallestTract;
  while(it != tracts.end())
  {
    TractType tract = *it;
    if(tract.size()<min)
    {
      smallestTract = tract;
      min = tract.size();
    }
    ++it;
  }






  it = tracts.begin();
  while(it != tracts.end())
  {
    TractType tract = *it;
    if(tract == smallestTract)
    {
      ++it;
      continue;
    }

    // Finding correspondences between points
    std::vector<int> correspondingIndices;
    TractType correspondingPoints;

    for(int i=0; i<smallestTract.size(); i++)
    {
      PointType p = smallestTract[i];

      double minDist = std::numeric_limits<float>::max();
      int correspondingIndex = 0;
      PointType correspondingPoint;

      // Search for the point on the second tract with the smallest distance
      // to p and memorize it
      for(int j=0; j<tract.size(); j++)
      {
        PointType p2 = tract[j];
        double dist = fabs( p.EuclideanDistanceTo(p2) );
        if(dist < minDist)
        {
          correspondingIndex = j;
          correspondingPoint = p2;
          minDist = dist;
        }
      }


      correspondingIndices.push_back(correspondingIndex);
      correspondingPoints.push_back(correspondingPoint);


    }

    if(smallestTract.size() != correspondingIndices.size())
    {

      MITK_ERROR << "smallest tract and correspondingIndices have no equal size";
      continue;
    }


    std::cout << "corresponding points\n";

    for(int i=0; i<smallestTract.size(); i++)
    {
      PointType p = smallestTract[i];
      PointType p2 = correspondingPoints[i];

      std::cout << "[" << p[0] << ", " << p[1] << ", " << p[2] << ", ["
                   << p2[0] << ", " << p2[1] << ", " << p2[2] << "]\n";
    }


    std::cout << std::endl;

    ++it;
  }




  // Get the values along the curves from a 3D images. should also contain info about the position on screen.
  std::vector< std::vector <double > > profiles;

  it = tracts.begin();
  while(it != tracts.end())
  {
    std::cout << "Tract\n";
    TractType tract = *it;
    TractType::iterator tractIt = tract.begin();

    std::vector<double> profile;

    while(tractIt != tract.end())
    {
      PointType p = *tractIt;
      std::cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';


      // Get value from image
      profile.push_back( (double)img->GetPixelValueByWorldCoordinate(p) );

      ++tractIt;
    }

    profiles.push_back(profile);
    std::cout << std::endl;

    ++it;
  }




  std::string title = "Fiber bundle plot";
  this->SetPlotTitle( title.c_str() );
  QPen pen( Qt::SolidLine );
  pen.setWidth(2);


  std::vector< std::vector<double> >::iterator profit = profiles.begin();

  int id=0;


  while(profit != profiles.end())
  {
    std::vector<double> profile = *profit;



    std::vector<double> xAxis;
    for(int i=0; i<profile.size(); ++i)
    {
      xAxis.push_back((double)i);
    }

    int curveId = this->InsertCurve( QString::number(id).toStdString().c_str() );
    this->SetCurveData( curveId, xAxis, profile );

    ++profit;
    id++;

  }

  this->Replot();










}


void QmitkTbssRoiAnalysisWidget::Boxplots()
{
  this->Clear();
}

void QmitkTbssRoiAnalysisWidget::drawBar(int x)
{

  m_Plot->detachItems(QwtPlotItem::Rtti_PlotMarker, true);


  QwtPlotMarker *mX = new QwtPlotMarker();
  //mX->setLabel(QString::fromLatin1("selected point"));
  mX->setLabelAlignment(Qt::AlignLeft | Qt::AlignBottom);
  mX->setLabelOrientation(Qt::Vertical);
  mX->setLineStyle(QwtPlotMarker::VLine);
  mX->setLinePen(QPen(Qt::black, 0, Qt::SolidLine));
  mX->setXValue(x);
  mX->attach(m_Plot);


  this->Replot();

}

QmitkTbssRoiAnalysisWidget::~QmitkTbssRoiAnalysisWidget()
{
  delete m_PlotPicker;

}




