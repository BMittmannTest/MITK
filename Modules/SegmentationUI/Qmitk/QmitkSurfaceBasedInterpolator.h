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

#ifndef QmitkSurfaceBasedInterpolator_h_Included
#define QmitkSurfaceBasedInterpolator_h_Included

#include "mitkSliceNavigationController.h"
#include "SegmentationUIExports.h"
#include "mitkDataNode.h"
#include "mitkDataStorage.h"
#include "mitkSurfaceInterpolationController.h"
#include "mitkToolManager.h"
//#include "mitkDiffSliceOperation.h"
//#include "mitkContourModelSet.h"
#include "mitkLabelSetImage.h"

#include <map>

#include <QWidget>

//For running 3D interpolation in background
#include <QtConcurrentRun>
#include <QFuture>
#include <QFutureWatcher>
#include <QTimer>

#include "ui_QmitkSurfaceBasedInterpolatorControls.h"

/**
  \brief GUI for surface-based interpolation.

  \ingroup ToolManagerEtAl
  \ingroup Widgets

  \sa QmitkInteractiveSegmentation
  \sa mitk::SurfaceInterpolationController

  There is a separate page describing the general design of QmitkInteractiveSegmentation: \ref QmitkInteractiveSegmentationTechnicalPage

  QmitkSurfaceBasedInterpolatorController is responsible to watch the GUI, to notice, which slice is currently
  visible. It triggers generation of interpolation suggestions and also triggers acception of suggestions.
*/

class SegmentationUI_EXPORT QmitkSurfaceBasedInterpolator : public QWidget
{
  Q_OBJECT

  public:

    QmitkSurfaceBasedInterpolator(QWidget* parent = 0, const char* name = 0);

    /**
      Initializes the widget. To be called once before real use.
    */
    void Initialize(mitk::DataStorage* storage);

    /**
      Removal of observers.
    */
    void Uninitialize();

    virtual ~QmitkSurfaceBasedInterpolator();

    void OnToolManagerWorkingDataModified();

   /**
      Just public because it is called by itk::Commands. You should not need to call this.
    */
    void OnSurfaceInterpolationInfoChanged(const itk::EventObject&);

    /**
     * @brief Set the visibility of the interpolation
     */
    void ShowInterpolationResult(bool);

  signals:

    void SignalShowMarkerNodes(bool);

  protected slots:

    void OnActivateWidget(bool);

    void OnAcceptInterpolationClicked();

    void OnSurfaceInterpolationFinished();

    void OnRunInterpolation();

    void OnShowMarkers(bool);

    void StartUpdateInterpolationTimer();

    void StopUpdateInterpolationTimer();

    void ChangeSurfaceColor();

private:

    void SetCurrentContourListID();

    mitk::SurfaceInterpolationController::Pointer m_SurfaceInterpolator;

    mitk::ToolManager::Pointer m_ToolManager;

    Ui::QmitkSurfaceBasedInterpolatorControls m_Controls;

    bool m_Initialized;

    bool m_Activated;

    unsigned int m_SurfaceInterpolationInfoChangedObserverTag;

//    mitk::DiffSliceOperation* m_doOperation;
//    mitk::DiffSliceOperation* m_undoOperation;

    mitk::DataNode::Pointer m_InterpolatedSurfaceNode;
    mitk::DataNode::Pointer m_3DContourNode;

    mitk::DataStorage::Pointer m_DataStorage;

    mitk::LabelSetImage::Pointer m_WorkingImage;

    QFuture<void> m_Future;
    QFutureWatcher<void> m_Watcher;
    QTimer* m_Timer;
};

#endif
