/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2010-03-31 16:40:27 +0200 (Mi, 31 Mrz 2010) $
Version:   $Revision: 21975 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/



#include <berryISelectionListener.h>
#include <berryIStructuredSelection.h>

#include <QmitkFunctionality.h>
#include "ui_QmitkDwiSoftwarePhantomViewControls.h"
#include <itkVectorImage.h>
#include <itkVectorContainer.h>
#include <itkOrientationDistributionFunction.h>

/*!
\brief QmitkDwiSoftwarePhantomView

\warning  This application module is not yet documented. Use "svn blame/praise/annotate" and ask the author to provide basic documentation.

\sa QmitkFunctionality
\ingroup Functionalities
*/

// Forward Qt class declarations


class QmitkDwiSoftwarePhantomView : public QmitkFunctionality
{

  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkDwiSoftwarePhantomView();
  virtual ~QmitkDwiSoftwarePhantomView();

  virtual void CreateQtPartControl(QWidget *parent);

  virtual void StdMultiWidgetAvailable (QmitkStdMultiWidget &stdMultiWidget);
  virtual void StdMultiWidgetNotAvailable();

  typedef itk::Image<unsigned char, 3>  ItkUcharImgType;
  typedef itk::Image<float, 3>          ItkFloatImgType;
  typedef itk::Vector<double,3>         GradientType;
  typedef std::vector<GradientType>     GradientListType;

  template<int ndirs> std::vector<itk::Vector<double,3> > MakeGradientList() ;

  protected slots:

  void GeneratePhantom();
  void OnSimulateBaselineToggle(int state);

protected:

  /// \brief called by QmitkFunctionality when DataManager's selection has changed
  virtual void OnSelectionChanged( std::vector<mitk::DataNode*> nodes );
  GradientListType GenerateHalfShell(int NPoints);

  Ui::QmitkDwiSoftwarePhantomViewControls* m_Controls;

  QmitkStdMultiWidget* m_MultiWidget;

  std::vector< mitk::DataNode::Pointer >    m_SignalRegionNodes;
  std::vector< ItkUcharImgType::Pointer >   m_SignalRegions;

  std::vector< QLabel* >    m_Labels;
  std::vector< QDoubleSpinBox* >    m_SpinFa;
  std::vector< QDoubleSpinBox* >    m_SpinAdc;
  std::vector< QDoubleSpinBox* >    m_SpinX;
  std::vector< QDoubleSpinBox* >    m_SpinY;
  std::vector< QDoubleSpinBox* >    m_SpinZ;
  std::vector< QDoubleSpinBox* >    m_SpinWeight;

  void UpdateGui();

private:

 };



