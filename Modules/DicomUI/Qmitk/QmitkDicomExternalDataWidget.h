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

#ifndef QmitkDicomExternalDataWidget_h
#define QmitkDicomExternalDataWidget_h

// #include <QmitkFunctionality.h>
#include "ui_QmitkDicomExternalDataWidgetControls.h"
#include "mitkDicomUIExports.h"

// include ctk
#include <ctkDICOMDatabase.h>
#include <ctkDICOMModel.h>
#include <ctkDICOMIndexer.h>
#include <ctkFileDialog.h>

//include QT
#include <QWidget>
#include <QString>
#include <QStringList>
#include <QModelIndex>
//For running dicom import in background
#include <QtConcurrentRun>
#include <QFuture>
#include <QFutureWatcher>
#include <QTimer>
#include <QProgressDialog>
#include <QLabel>

/*!
\brief QmitkDicomExternalDataWidget 

\warning  This application module is not yet documented. Use "svn blame/praise/annotate" and ask the author to provide basic documentation.

\sa QmitkFunctionality
\ingroup Functionalities
*/
class MITK_DICOMUI_EXPORT QmitkDicomExternalDataWidget : public QWidget
{  
   // this is needed for all Qt objects that should have a Qt meta-object
   // (everything that derives from QObject and wants to have signal/slots)
   Q_OBJECT

public:  

   static const std::string Widget_ID;

   QmitkDicomExternalDataWidget(QWidget *parent);
   virtual ~QmitkDicomExternalDataWidget();

   virtual void CreateQtPartControl(QWidget *parent);

   /* @brief   Initializes the widget. This method has to be called before widget can start. 
   * @param dataStorage The data storage the widget should work with.
   * @param multiWidget The corresponding multiwidget were the ct Image is displayed and the user should define his path.
   * @param imageNode  The image node which will be the base of mitral processing
   */
   void Initialize();

signals:
    void SignalChangePage(int);
    void SignalAddDicomData(const QString&);
    void SignalAddDicomData(const QStringList&);
    void SignalDicomToDataManager(const QStringList&);

   public slots:

       /// @brief Called when import CD or import Folder was clicked.
       void OnFolderCDImport();

       /// @brief Called when import directory was selected.
       void OnFileSelectedAddExternalData(QString);

       /// @brief Called when download button was clicked.
       void OnDownloadButtonClicked();

       /// @brief Called when view button was clicked.
       void OnViewButtonClicked();
    
       /// @brief Called when cancel button was clicked.
       void OnCancelButtonClicked();

       /// @brief   Called when search parameters change.
        void OnSearchParameterChanged();

        /// @brief   Called when import progress change.
        void OnProgress(int progress);

protected:

    /// \brief Get the list of filepath from current selected index in TreeView. All file paths referring to the index will be returned.
    void GetFileNamesFromIndex(QStringList& filePaths);

    void AddDicomTemporary(QString directory);

    ctkDICOMDatabase* m_ExternalDatabase;
    ctkDICOMModel* m_ExternalModel;
    ctkDICOMIndexer* m_ExternalIndexer;
    
    ctkFileDialog* m_ImportDialog;
    QProgressDialog* m_ProgressDialog;
    QLabel* m_ProgressDialogLabel;

    Ui::QmitkDicomExternalDataWidgetControls* m_Controls;

    QFuture<void> m_Future;
    QFutureWatcher<void> m_Watcher;
    QTimer* m_Timer;
    QString* m_DirectoryName;

};



#endif // _QmitkDicomExternalDataWidget_H_INCLUDED

