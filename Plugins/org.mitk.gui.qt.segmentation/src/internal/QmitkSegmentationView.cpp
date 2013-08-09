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

#include "mitkProperties.h"
#include "mitkSegTool2D.h"

#include "QmitkStdMultiWidget.h"
#include "QmitkNewSegmentationDialog.h"

#include <QMessageBox>
#include <QCompleter>
#include <QStringListModel>
#include <QInputDialog>
#include <QFileDialog>
#include <QDateTime>

#include <berryIWorkbenchPage.h>

#include "QmitkSegmentationView.h"
#include "QmitkSegmentationOrganNamesHandling.cpp"


#include <mitkSurfaceToImageFilter.h>

#include "mitkVtkResliceInterpolationProperty.h"
#include "mitkRenderingModeProperty.h"
#include "mitkGetModuleContext.h"
#include "mitkModule.h"
#include "mitkModuleRegistry.h"

#include "mitkLabelSetImage.h"
#include "mitkColormapProperty.h"
#include "mitkImageWriteAccessor.h"
#include "mitkNrrdLabelSetImageWriter.h"
#include "mitkNrrdLabelSetImageReader.h"

#include "mitkModuleResource.h"
#include "mitkStatusBar.h"
#include "mitkApplicationCursor.h"

#include "mitkSegmentationObjectFactory.h"

const std::string QmitkSegmentationView::VIEW_ID = "org.mitk.views.segmentation";

// public methods

QmitkSegmentationView::QmitkSegmentationView()
:m_Parent(NULL)
,m_Controls(NULL)
,m_MultiWidget(NULL)
,m_DataSelectionChanged(false)
,m_MouseCursorSet(false)
{
  RegisterSegmentationObjectFactory();
  mitk::NodePredicateDataType::Pointer isDwi = mitk::NodePredicateDataType::New("DiffusionImage");
  mitk::NodePredicateDataType::Pointer isDti = mitk::NodePredicateDataType::New("TensorImage");
  mitk::NodePredicateDataType::Pointer isQbi = mitk::NodePredicateDataType::New("QBallImage");

  mitk::NodePredicateOr::Pointer isDiffusionImage = mitk::NodePredicateOr::New();
  isDiffusionImage->AddPredicate(isDwi);
  isDiffusionImage->AddPredicate(isDti);
  isDiffusionImage->AddPredicate(isQbi);

  m_IsOfTypeImagePredicate = mitk::NodePredicateOr::New(isDiffusionImage, mitk::TNodePredicateDataType<mitk::Image>::New());

  m_IsOfTypeLabelSetImagePredicate = mitk::NodePredicateDataType::New("LabelSetImage");
  m_IsNotLabelSetImagePredicate = mitk::NodePredicateAnd::New();
  m_IsNotLabelSetImagePredicate->AddPredicate( m_IsOfTypeImagePredicate );
  m_IsNotLabelSetImagePredicate->AddPredicate( mitk::NodePredicateNot::New( m_IsOfTypeLabelSetImagePredicate));
}

QmitkSegmentationView::~QmitkSegmentationView()
{
  delete m_Controls;
}

void QmitkSegmentationView::NewNodeObjectsGenerated(mitk::ToolManager::DataVectorType* nodes)
{
  if (!nodes) return;

  mitk::ToolManager* toolManager = m_Controls->m_ManualToolSelectionBox2D->GetToolManager();
  if (!toolManager) return;
  for (mitk::ToolManager::DataVectorType::iterator iter = nodes->begin(); iter != nodes->end(); ++iter)
  {
    this->FireNodeSelected( *iter );
    // only last iteration meaningful, multiple generated objects are not taken into account here
  }
}

void QmitkSegmentationView::Visible()
{
  if (m_DataSelectionChanged)
  {
    this->OnSelectionChanged(this->GetDataManagerSelection());
  }
}

void QmitkSegmentationView::Activated()
{
  // should be moved to ::BecomesVisible() or similar
  if( m_Controls )
  {
    //m_Controls->m_ManualToolSelectionBox2D->SetAutoShowNamesWidth(m_Controls->m_ManualToolSelectionBox2D->minimumSizeHint().width()+1);
    m_Controls->m_ManualToolSelectionBox2D->SetAutoShowNamesWidth(250);
    m_Controls->m_ManualToolSelectionBox2D->setEnabled( true );
    //m_Controls->m_ManualToolSelectionBox3D->SetAutoShowNamesWidth(m_Controls->m_ManualToolSelectionBox3D->minimumSizeHint().width()+1);
    m_Controls->m_ManualToolSelectionBox3D->SetAutoShowNamesWidth(260);
    m_Controls->m_ManualToolSelectionBox3D->setEnabled( true );
//    m_Controls->m_OrganToolSelectionBox->setEnabled( true );
//    m_Controls->m_LesionToolSelectionBox->setEnabled( true );

//    m_Controls->m_SlicesInterpolator->Enable3DInterpolation( m_Controls->widgetStack->currentWidget() == m_Controls->pageManual );

    mitk::DataStorage::SetOfObjects::ConstPointer segmentations = this->GetDefaultDataStorage()->GetSubset( m_IsOfTypeLabelSetImagePredicate );

    mitk::DataStorage::SetOfObjects::ConstPointer image = this->GetDefaultDataStorage()->GetSubset( m_IsNotLabelSetImagePredicate );
    if (!image->empty()) {
      OnSelectionChanged(*image->begin());
    }

  }
  this->SetToolManagerSelection(m_Controls->patImageSelector->GetSelectedNode(), m_Controls->segImageSelector->GetSelectedNode());
}

void QmitkSegmentationView::Deactivated()
{
  if( m_Controls )
  {
    m_Controls->m_ManualToolSelectionBox2D->setEnabled( false );
    m_Controls->m_ManualToolSelectionBox3D->setEnabled( false );
    //deactivate all tools
    m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->ActivateTool(-1);
    m_Controls->m_SlicesInterpolator->EnableInterpolation( false );

    // gets the context of the "Mitk" (Core) module (always has id 1)
    // TODO Workaround until CTK plugincontext is available
    mitk::ModuleContext* context = mitk::ModuleRegistry::GetModule(1)->GetModuleContext();
    // Workaround end
    mitk::ServiceReference serviceRef = context->GetServiceReference<mitk::PlanePositionManagerService>();
    mitk::PlanePositionManagerService* service = dynamic_cast<mitk::PlanePositionManagerService*>(context->GetService(serviceRef));
    service->RemoveAllPlanePositions();
  }
}

void QmitkSegmentationView::StdMultiWidgetAvailable( QmitkStdMultiWidget& stdMultiWidget )
{
  SetMultiWidget(&stdMultiWidget);
}

void QmitkSegmentationView::StdMultiWidgetNotAvailable()
{
  SetMultiWidget(NULL);
}

void QmitkSegmentationView::StdMultiWidgetClosed( QmitkStdMultiWidget& /*stdMultiWidget*/ )
{
  SetMultiWidget(NULL);
}

void QmitkSegmentationView::SetMultiWidget(QmitkStdMultiWidget* multiWidget)
{
  // save the current multiwidget as the working widget
  m_MultiWidget = multiWidget;

  if (m_Parent)
  {
    m_Parent->setEnabled(m_MultiWidget);
  }

  // tell the interpolation about toolmanager and multiwidget (and data storage)
  if (m_Controls && m_MultiWidget)
  {
    mitk::ToolManager* toolManager = m_Controls->m_ManualToolSelectionBox2D->GetToolManager();
    m_Controls->m_SlicesInterpolator->SetDataStorage( this->GetDefaultDataStorage());
    QList<mitk::SliceNavigationController*> controllers;
    controllers.push_back(m_MultiWidget->GetRenderWindow1()->GetSliceNavigationController());
    controllers.push_back(m_MultiWidget->GetRenderWindow2()->GetSliceNavigationController());
    controllers.push_back(m_MultiWidget->GetRenderWindow3()->GetSliceNavigationController());
    m_Controls->m_SlicesInterpolator->Initialize( toolManager, controllers );
  }
}

void QmitkSegmentationView::OnPreferencesChanged(const berry::IBerryPreferences* prefs)
{
  m_AutoSelectionEnabled = prefs->GetBool("auto selection", false);
}

void QmitkSegmentationView::OnLabelListModified(const QStringList& list)
{
   QStringListModel* completeModel = static_cast<QStringListModel*> (m_Completer->model());
   completeModel->setStringList(list);
}

void QmitkSegmentationView::OnNewLabel()
{
    // ask about the name and organ type of the new segmentation
    QmitkNewSegmentationDialog* dialog = new QmitkNewSegmentationDialog( m_Parent );

    QString storedList = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
    QStringList organColors;
    if (storedList.isEmpty())
    {
        organColors = GetDefaultOrganColorString();
    }
    else
    {
        /*
        a couple of examples of how organ names are stored:

        a simple item is built up like 'name#AABBCC' where #AABBCC is the hexadecimal notation of a color as known from HTML

        items are stored separated by ';'
        this makes it necessary to escape occurrences of ';' in name.
        otherwise the string "hugo;ypsilon#AABBCC;eugen#AABBCC" could not be parsed as two organs
        but we would get "hugo" and "ypsilon#AABBCC" and "eugen#AABBCC"

        so the organ name "hugo;ypsilon" is stored as "hugo\;ypsilon"
        and must be unescaped after loading

        the following lines could be one split with Perl's negative lookbehind
        */

        // recover string list from BlueBerry view's preferences
        QString storedString = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
        MITK_DEBUG << "storedString: " << storedString.toStdString();
        // match a string consisting of any number of repetitions of either "anything but ;" or "\;". This matches everything until the next unescaped ';'
        QRegExp onePart("(?:[^;]|\\\\;)*");
        MITK_DEBUG << "matching " << onePart.pattern().toStdString();
        int count = 0;
        int pos = 0;
        while( (pos = onePart.indexIn( storedString, pos )) != -1 )
        {
            ++count;
            int length = onePart.matchedLength();
            if (length == 0) break;
            QString matchedString = storedString.mid(pos, length);
            MITK_DEBUG << "   Captured length " << length << ": " << matchedString.toStdString();
            pos += length + 1; // skip separating ';'

            // unescape possible occurrences of '\;' in the string
            matchedString.replace("\\;", ";");

            // add matched string part to output list
            organColors << matchedString;
        }
        MITK_DEBUG << "Captured " << count << " organ name/colors";
    }

    dialog->SetSuggestionList( organColors );

    int dialogReturnValue = dialog->exec();

    if ( dialogReturnValue == QDialog::Rejected ) return; // user clicked cancel or pressed Esc or something similar

    mitk::Color color = dialog->GetColor();

    mitk::DataNode* workingNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetWorkingData(0);
    if (!workingNode) return;

    mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>(workingNode->GetData());
    if ( !lsImage ) return;

    lsImage->AddLabel(dialog->GetSegmentationName().toStdString(), color);
}

void QmitkSegmentationView::OnCreateSurface(int index)
{
}

void QmitkSegmentationView::OnCombineAndCreateSurface( const QList<QTableWidgetSelectionRange>& ranges )
{

}

void QmitkSegmentationView::OnCombineAndCreateMask( const QList<QTableWidgetSelectionRange>& ranges )
{

}

void QmitkSegmentationView::OnRenameLabel(int index)
{
    // ask about the name and organ type of the new segmentation
    QmitkNewSegmentationDialog* dialog = new QmitkNewSegmentationDialog( m_Parent ); // needs a QWidget as parent, "this" is not QWidget

    QString storedList = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
    QStringList organColors;
    if (storedList.isEmpty())
    {
        organColors = GetDefaultOrganColorString();
    }
    else
    {
        /*
        a couple of examples of how organ names are stored:

        a simple item is built up like 'name#AABBCC' where #AABBCC is the hexadecimal notation of a color as known from HTML

        items are stored separated by ';'
        this makes it necessary to escape occurrences of ';' in name.
        otherwise the string "hugo;ypsilon#AABBCC;eugen#AABBCC" could not be parsed as two organs
        but we would get "hugo" and "ypsilon#AABBCC" and "eugen#AABBCC"

        so the organ name "hugo;ypsilon" is stored as "hugo\;ypsilon"
        and must be unescaped after loading

        the following lines could be one split with Perl's negative lookbehind
        */

        // recover string list from BlueBerry view's preferences
        QString storedString = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
        MITK_DEBUG << "storedString: " << storedString.toStdString();
        // match a string consisting of any number of repetitions of either "anything but ;" or "\;". This matches everything until the next unescaped ';'
        QRegExp onePart("(?:[^;]|\\\\;)*");
        MITK_DEBUG << "matching " << onePart.pattern().toStdString();
        int count = 0;
        int pos = 0;
        while( (pos = onePart.indexIn( storedString, pos )) != -1 )
        {
            ++count;
            int length = onePart.matchedLength();
            if (length == 0) break;
            QString matchedString = storedString.mid(pos, length);
            MITK_DEBUG << "   Captured length " << length << ": " << matchedString.toStdString();
            pos += length + 1; // skip separating ';'

            // unescape possible occurrences of '\;' in the string
            matchedString.replace("\\;", ";");

            // add matched string part to output list
            organColors << matchedString;
        }
        MITK_DEBUG << "Captured " << count << " organ name/colors";
    }

    dialog->SetSuggestionList( organColors );

    int dialogReturnValue = dialog->exec();

    if ( dialogReturnValue == QDialog::Rejected ) return; // user clicked cancel or pressed Esc or something similar

    mitk::DataNode* workingNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetWorkingData(0);
    if (!workingNode) return;

    mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>(workingNode->GetData());
    if ( !lsImage ) return;

    lsImage->RenameLabel(index, dialog->GetSegmentationName().toStdString(), dialog->GetColor());
}

void QmitkSegmentationView::OnWorkingNodeVisibilityChanged()
{
  m_Controls->m_ManualToolSelectionBox2D->setEnabled(false);
  if (!m_Controls->m_ManualToolSelectionBox2D->isEnabled())
  {
    this->UpdateWarningLabel("The selected segmentation is currently not visible!");
  }
  else
  {
    this->UpdateWarningLabel("");
  }
}

void QmitkSegmentationView::NodeAdded(const mitk::DataNode *node)
{

}

void QmitkSegmentationView::NodeRemoved(const mitk::DataNode* node)
{
  bool isHelperObject(false);
  node->GetBoolProperty("helper object", isHelperObject);

  mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>(node->GetData());

  if(lsImage && !isHelperObject)
  {
    //First of all remove all possible contour markers of the segmentation
    mitk::DataStorage::SetOfObjects::ConstPointer allContourMarkers = this->GetDataStorage()->GetDerivations(node, mitk::NodePredicateProperty::New("isContourMarker"
                                                                            , mitk::BoolProperty::New(true)));

    // gets the context of the "Mitk" (Core) module (always has id 1)
    // TODO Workaround until CTK plugincontext is available
    mitk::ModuleContext* context = mitk::ModuleRegistry::GetModule(1)->GetModuleContext();
    // Workaround end
    mitk::ServiceReference serviceRef = context->GetServiceReference<mitk::PlanePositionManagerService>();

    mitk::PlanePositionManagerService* service = dynamic_cast<mitk::PlanePositionManagerService*>(context->GetService(serviceRef));

    for (mitk::DataStorage::SetOfObjects::ConstIterator it = allContourMarkers->Begin(); it != allContourMarkers->End(); ++it)
    {
      std::string nodeName = node->GetName();
      unsigned int t = nodeName.find_last_of(" ");
      unsigned int id = atof(nodeName.substr(t+1).c_str())-1;

      service->RemovePlanePosition(id);

      this->GetDataStorage()->Remove(it->Value());
    }

    if ((m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetWorkingData(0) == node) && m_Controls->patImageSelector->GetSelectedNode().IsNotNull())
    {
      this->SetToolManagerSelection(m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetReferenceData(0), NULL);
      this->UpdateWarningLabel("Create a segmentation!");
    }
    mitk::SurfaceInterpolationController::GetInstance()->RemoveSegmentationFromContourList(lsImage);
  }

  if((m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetReferenceData(0) == node))
  {
    //as we don't know which node was actually removed e.g. our reference node, disable 'New Segmentation' button.
    //consider the case that there is no more image in the datastorage
    this->SetToolManagerSelection(NULL, NULL);
  }
}

void QmitkSegmentationView::OnSearchLabel()
{
    m_Controls->m_LabelSetTableWidget->SetActiveLabel(m_Controls->m_LabelSearchBox->text().toStdString());
}

void QmitkSegmentationView::OnLoadLabelSet()
{
    mitk::DataNode* refNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetReferenceData(0);
    if (!refNode) return;

    mitk::Image* refImage = dynamic_cast<mitk::Image*>( refNode->GetData() );
    if (!refImage) return;

    std::string fileExtensions("Segmentation files (*.lset);;");
    QString qfileName = QFileDialog::getOpenFileName(m_Parent, "Open segmentation", "", fileExtensions.c_str() );
    if (qfileName.isEmpty() ) return;

    mitk::NrrdLabelSetImageReader<unsigned char>::Pointer reader = mitk::NrrdLabelSetImageReader<unsigned char>::New();
    reader->SetFileName(qfileName.toLatin1());

    this->WaitCursorOn();

    try
    {
        reader->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
        MITK_ERROR << "Exception caught: " << excep.GetDescription();
        QMessageBox::information(m_Parent, "Load segmentation", "Could not load active segmentation. See error log for details.\n");
        this->WaitCursorOff();
    }

    this->WaitCursorOff();

    mitk::LabelSetImage::Pointer lsImage = reader->GetOutput();

    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    newNode->SetData(lsImage);
    newNode->SetName( lsImage->GetLabelSetName() );

    this->GetDataStorage()->Add(newNode, refNode);

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkSegmentationView::OnSaveLabelSet()
{
  mitk::DataNode* workingNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetWorkingData(0);
  if (!workingNode) return;

  mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>(workingNode->GetData());
  if ( !lsImage ) return;

  // update the name in case has changed
  lsImage->SetLabelSetName( workingNode->GetName() );

  //set last modification date and time
  QDateTime cTime = QDateTime::currentDateTime ();
  lsImage->SetLabelSetLastModified( cTime.toString().toStdString() );

  QString proposedFileName = QString::fromStdString(workingNode->GetName());
  QString selected_suffix("segmentation files (*.lset)");
  QString qfileName = QFileDialog::getSaveFileName(m_Parent, QString("Save segmentation..."), proposedFileName, selected_suffix);
  if (qfileName.isEmpty() ) return;

  mitk::NrrdLabelSetImageWriter<unsigned char>::Pointer writer = mitk::NrrdLabelSetImageWriter<unsigned char>::New();
  writer->SetFileName(qfileName.toStdString());
  std::string newName = itksys::SystemTools::GetFilenameWithoutExtension(qfileName.toStdString());
  lsImage->SetLabelSetName(newName);
  workingNode->SetName(newName);
  writer->SetInput(lsImage);

  this->WaitCursorOn();

  try
  {
      writer->Update();
  }
  catch ( itk::ExceptionObject & excep )
  {
    MITK_ERROR << "Exception caught: " << excep.GetDescription();
    QMessageBox::information(m_Parent, "Save segmenation", "Could not save active segmentation. See error log for details.\n\n");
    this->WaitCursorOff();
  }

  this->WaitCursorOff();
}

void QmitkSegmentationView::OnNewLabelSet()
{
  mitk::DataNode::Pointer refNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetReferenceData(0);
  if (refNode.IsNotNull())
  {
    mitk::Image::Pointer refImage = dynamic_cast<mitk::Image*>( refNode->GetData() );
    if (refImage.IsNotNull())
    {
        bool ok = false;
        QString refNodeName = QString::fromStdString(refNode->GetName());
        refNodeName.append("-");
        QString newName = QInputDialog::getText(m_Parent, "New segmentation", "Set name:", QLineEdit::Normal, refNodeName, &ok);
        if (!ok) return;

        mitk::PixelType pixelType(mitk::MakeScalarPixelType<unsigned char>() );
        mitk::LabelSetImage::Pointer labelSetImage = mitk::LabelSetImage::New();

        if (refImage->GetDimension() == 2)
        {
            const unsigned int dimensions[] = { refImage->GetDimension(0), refImage->GetDimension(1), 1 };
            labelSetImage->Initialize(pixelType, 3, dimensions);
        }
        else
        {
            labelSetImage->Initialize(pixelType, refImage->GetDimension(), refImage->GetDimensions());
        }

        unsigned int byteSize = sizeof(unsigned char);
        for (unsigned int dim = 0; dim < labelSetImage->GetDimension(); ++dim)
        {
            byteSize *= labelSetImage->GetDimension(dim);
        }

        mitk::ImageWriteAccessor accessor(static_cast<mitk::Image*>(labelSetImage));
        memset( accessor.GetData(), 0, byteSize );

        mitk::TimeSlicedGeometry::Pointer originalGeometry = refImage->GetTimeSlicedGeometry()->Clone();
        labelSetImage->SetGeometry( originalGeometry );

        mitk::DataNode::Pointer newNode = mitk::DataNode::New();
        newNode->SetData(labelSetImage);
        newNode->SetName(newName.toStdString());
        labelSetImage->SetLabelSetName(newName.toStdString());

        this->GetDataStorage()->Add(newNode, refNode);

        mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    }
  }
}

void QmitkSegmentationView::CreateLabelFromSurface()
{
    /*
  mitk::DataNode* imageNode
      = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetReferenceData(0);

  if (!imageNode) return;

  mitk::DataNode* surfaceNode = m_Controls->MaskSurfaces->GetSelectedNode();
  if (!surfaceNode) return;

  mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>( surfaceNode->GetData() );
  if (surface.IsNull()) return;

//adapt this one or write a new one for label set images
  mitk::SurfaceToImageFilter::Pointer s2iFilter
      = mitk::SurfaceToImageFilter::New();

  s2iFilter->MakeOutputBinaryOn();
  s2iFilter->SetInput(surface);
  s2iFilter->SetImage(image);
  s2iFilter->Update();
*/
// ask about the name and organ type of the new segmentation
    QmitkNewSegmentationDialog* dialog = new QmitkNewSegmentationDialog( m_Parent );

    QString storedList = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
    QStringList organColors;
    if (storedList.isEmpty())
    {
        organColors = GetDefaultOrganColorString();
    }
    else
    {
        // recover string list from BlueBerry view's preferences
        QString storedString = QString::fromStdString( this->GetPreferences()->GetByteArray("Organ-Color-List","") );
        // match a string consisting of any number of repetitions of either "anything but ;" or "\;". This matches everything until the next unescaped ';'
        QRegExp onePart("(?:[^;]|\\\\;)*");
        int count = 0;
        int pos = 0;
        while( (pos = onePart.indexIn( storedString, pos )) != -1 )
        {
            ++count;
            int length = onePart.matchedLength();
            if (length == 0) break;
            QString matchedString = storedString.mid(pos, length);
            MITK_DEBUG << "   Captured length " << length << ": " << matchedString.toStdString();
            pos += length + 1; // skip separating ';'

            // unescape possible occurrences of '\;' in the string
            matchedString.replace("\\;", ";");

            // add matched string part to output list
            organColors << matchedString;
        }
    }

    dialog->SetSuggestionList( organColors );

    int dialogReturnValue = dialog->exec();

    if ( dialogReturnValue == QDialog::Rejected ) return;

    mitk::Color color = dialog->GetColor();

    mitk::DataNode* workingNode = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetWorkingData(0);
    if (!workingNode) return;

    mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>(workingNode->GetData());
    if ( !lsImage ) return;

    lsImage->AddLabel(dialog->GetSegmentationName().toStdString(), color);
}

// protected

void QmitkSegmentationView::OnPatientComboBoxSelectionChanged( const mitk::DataNode* node )
{
  //mitk::DataNode* selectedNode = const_cast<mitk::DataNode*>(node);
  if( node != NULL )
  {
    m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->ActivateTool(-1);
    this->UpdateWarningLabel("");
    mitk::DataNode* segNode = m_Controls->segImageSelector->GetSelectedNode();
    if (segNode)
    {
      mitk::DataStorage::SetOfObjects::ConstPointer possibleParents = this->GetDefaultDataStorage()->GetSources( segNode, m_IsNotLabelSetImagePredicate );
      bool isSourceNode(false);

      for (mitk::DataStorage::SetOfObjects::ConstIterator it = possibleParents->Begin(); it != possibleParents->End(); it++)
      {
        if (it.Value() == node)
          isSourceNode = true;
      }

      if ( !isSourceNode && (!this->CheckForSameGeometry(segNode, node) || possibleParents->Size() > 0 ))
      {
        this->SetToolManagerSelection(const_cast<mitk::DataNode*>(node), NULL);
        this->UpdateWarningLabel("The selected patient image does not\n match with the selected segmentation!");
      }
      else if ((!isSourceNode && this->CheckForSameGeometry(segNode, node)) || isSourceNode )
      {
        this->SetToolManagerSelection(const_cast<mitk::DataNode*>(node), segNode);
        //Doing this we can assure that the segmenation is always visible if the segmentation and the patient image are
        //loaded separately
        int layer(10);
        node->GetIntProperty("layer", layer);
        layer++;
        segNode->SetProperty("layer", mitk::IntProperty::New(layer));
        this->UpdateWarningLabel("");
      }
    }
    else
    {
      this->SetToolManagerSelection(const_cast<mitk::DataNode*>(node), NULL);
      this->UpdateWarningLabel("Create a segmentation");
    }
  }
  else
  {
    this->UpdateWarningLabel("Please load an image!");
  }
}

void QmitkSegmentationView::OnSegmentationComboBoxSelectionChanged(const mitk::DataNode *node)
{
  if ( !node)
    return;

  mitk::DataNode* refNode = m_Controls->patImageSelector->GetSelectedNode();
  mitk::DataNode* segNode = const_cast<mitk::DataNode*>(node);

  m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->ActivateTool(-1);

  if (m_AutoSelectionEnabled)
  {
    this->OnSelectionChanged(segNode);
  }
  else
  {
    mitk::DataStorage::SetOfObjects::ConstPointer possibleParents = this->GetDefaultDataStorage()->GetSources( segNode, m_IsNotLabelSetImagePredicate );

    if ( possibleParents->Size() == 1 )
    {
      mitk::DataNode* parentNode = possibleParents->ElementAt(0);

      if (parentNode != refNode)
      {
        this->UpdateWarningLabel("The selected segmentation does not\n match with the selected patient image!");
        this->SetToolManagerSelection(NULL, segNode);
      }
      else
      {
        this->UpdateWarningLabel("");
        segNode->SetVisibility(true);
        this->SetToolManagerSelection(refNode, segNode);
        mitk::DataStorage::SetOfObjects::ConstPointer otherSegmentations = this->GetDataStorage()->GetSubset(m_IsOfTypeLabelSetImagePredicate);
        for(mitk::DataStorage::SetOfObjects::const_iterator iter = otherSegmentations->begin(); iter != otherSegmentations->end(); ++iter)
        {
          mitk::DataNode* _otherSegNode = *iter;
          if (dynamic_cast<mitk::Image*>(_otherSegNode->GetData()) != segNode->GetData())
            _otherSegNode->SetVisibility(false);
        }
        mitk::RenderingManager::GetInstance()->RequestUpdateAll();
      }
    }
    /* // todo: ask for a use case
    else if (refNode && this->CheckForSameGeometry(segNode, refNode))
    {
      this->UpdateWarningLabel("");
      this->SetToolManagerSelection(refNode, segNode);
    }
    */
    else if (!refNode || !this->CheckForSameGeometry(node, refNode))
    {
      this->UpdateWarningLabel("Please select or load the right patient image!");
    }
  }
  if (!node->IsVisible(mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget1"))))
    this->UpdateWarningLabel("The selected segmentation is currently not visible!");
}

void QmitkSegmentationView::OnShowMarkerNodes (bool state)
{
  mitk::SegTool2D::Pointer manualSegmentationTool;

  unsigned int numberOfExistingTools = m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetTools().size();

  for(unsigned int i = 0; i < numberOfExistingTools; i++)
  {
    manualSegmentationTool = dynamic_cast<mitk::SegTool2D*>(m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->GetToolById(i));

    if (manualSegmentationTool)
    {
      if(state == true)
      {
        manualSegmentationTool->SetShowMarkerNodes( true );
      }
      else
      {
        manualSegmentationTool->SetShowMarkerNodes( false );
      }
    }
  }
}


void QmitkSegmentationView::OnSelectionChanged(mitk::DataNode* node)
{
  std::vector<mitk::DataNode*> nodes;
  nodes.push_back( node );
  this->OnSelectionChanged( nodes );
}

void QmitkSegmentationView::OnSelectionChanged(std::vector<mitk::DataNode*> nodes)
{

}

void QmitkSegmentationView::OnContourMarkerSelected(const mitk::DataNode *node)
{
  QmitkRenderWindow* selectedRenderWindow = 0;
  QmitkRenderWindow* RenderWindow1 =
      this->GetActiveStdMultiWidget()->GetRenderWindow1();
  QmitkRenderWindow* RenderWindow2 =
      this->GetActiveStdMultiWidget()->GetRenderWindow2();
  QmitkRenderWindow* RenderWindow3 =
      this->GetActiveStdMultiWidget()->GetRenderWindow3();
  QmitkRenderWindow* RenderWindow4 =
      this->GetActiveStdMultiWidget()->GetRenderWindow4();
  bool PlanarFigureInitializedWindow = false;

  // find initialized renderwindow
  if (node->GetBoolProperty("PlanarFigureInitializedWindow",
                PlanarFigureInitializedWindow, RenderWindow1->GetRenderer()))
  {
    selectedRenderWindow = RenderWindow1;
  }
  if (!selectedRenderWindow && node->GetBoolProperty(
        "PlanarFigureInitializedWindow", PlanarFigureInitializedWindow,
        RenderWindow2->GetRenderer()))
  {
    selectedRenderWindow = RenderWindow2;
  }
  if (!selectedRenderWindow && node->GetBoolProperty(
        "PlanarFigureInitializedWindow", PlanarFigureInitializedWindow,
        RenderWindow3->GetRenderer()))
  {
    selectedRenderWindow = RenderWindow3;
  }
  if (!selectedRenderWindow && node->GetBoolProperty(
        "PlanarFigureInitializedWindow", PlanarFigureInitializedWindow,
        RenderWindow4->GetRenderer()))
  {
    selectedRenderWindow = RenderWindow4;
  }

  // make node visible
  if (selectedRenderWindow)
  {
    std::string nodeName = node->GetName();
    unsigned int t = nodeName.find_last_of(" ");
    unsigned int id = atof(nodeName.substr(t+1).c_str())-1;

    // gets the context of the "Mitk" (Core) module (always has id 1)
    // TODO Workaround until CTL plugincontext is available
    mitk::ModuleContext* context = mitk::ModuleRegistry::GetModule(1)->GetModuleContext();
    // Workaround end
    mitk::ServiceReference serviceRef = context->GetServiceReference<mitk::PlanePositionManagerService>();

    mitk::PlanePositionManagerService* service = dynamic_cast<mitk::PlanePositionManagerService*>(context->GetService(serviceRef));
    selectedRenderWindow->GetSliceNavigationController()->ExecuteOperation(service->GetPlanePosition(id));
    selectedRenderWindow->GetRenderer()->GetDisplayGeometry()->Fit();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}

void QmitkSegmentationView::OnTabWidgetChanged(int id)
{
  //2D Tab ID = 0
  //3D Tab ID = 1
  if (id == 0)
  {
    //Hide 3D selection box, show 2D selection box
    m_Controls->m_ManualToolSelectionBox3D->hide();
    m_Controls->m_ManualToolSelectionBox2D->show();
    //Deactivate possible active tool
    m_Controls->m_ManualToolSelectionBox3D->GetToolManager()->ActivateTool(-1);

    //TODO Remove possible visible interpolations -> Maybe changes in SlicesInterpolator
  }
  else
  {
    //Hide 3D selection box, show 2D selection box
    m_Controls->m_ManualToolSelectionBox2D->hide();
    m_Controls->m_ManualToolSelectionBox3D->show();
    //Deactivate possible active tool
    m_Controls->m_ManualToolSelectionBox2D->GetToolManager()->ActivateTool(-1);
  }
}


void QmitkSegmentationView::SetToolManagerSelection(mitk::DataNode* referenceData, mitk::DataNode* workingData)
{
  // called as a result of new BlueBerry selections
  //   tells the ToolManager for manual segmentation about new selections
  //   updates GUI information about what the user should select
  mitk::ToolManager* toolManager = m_Controls->m_ManualToolSelectionBox2D->GetToolManager();
  toolManager->SetReferenceData(referenceData);
  toolManager->SetWorkingData(workingData);

  if (workingData)
  {
      mitk::LabelSetImage* lsImage = dynamic_cast<mitk::LabelSetImage*>( workingData->GetData());
      m_Controls->m_LabelSetTableWidget->SetActiveLabelSetImage(lsImage);
  }
  else
      m_Controls->m_LabelSetTableWidget->SetActiveLabelSetImage(NULL);

  // check original image
  m_Controls->m_btNewLabelSet->setEnabled(referenceData != NULL);
  m_Controls->m_btLoadLabelSet->setEnabled(referenceData != NULL);
  if (referenceData)
  {
    this->UpdateWarningLabel("");
    disconnect( m_Controls->patImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
          this, SLOT( OnPatientComboBoxSelectionChanged( const mitk::DataNode* ) ) );

    m_Controls->patImageSelector->setCurrentIndex( m_Controls->patImageSelector->Find(referenceData) );

    connect( m_Controls->patImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
         this, SLOT( OnPatientComboBoxSelectionChanged( const mitk::DataNode* ) ) );
  }

  // check segmentation
  if (referenceData)
  {
    if (workingData)
    {
      this->FireNodeSelected(const_cast<mitk::DataNode*>(workingData));
//      mitk::RenderingManager::GetInstance()->InitializeViews(workingData->GetData()->GetTimeSlicedGeometry(),
//                                                             mitk::RenderingManager::REQUEST_UPDATE_ALL, true );

//      if( m_Controls->widgetStack->currentIndex() == 0 )
//      {
        disconnect( m_Controls->segImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
              this, SLOT( OnSegmentationComboBoxSelectionChanged( const mitk::DataNode* ) ) );

        m_Controls->segImageSelector->setCurrentIndex(m_Controls->segImageSelector->Find(workingData));

        connect( m_Controls->segImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
             this, SLOT( OnSegmentationComboBoxSelectionChanged(const mitk::DataNode*)) );
//      }
    }
    m_Controls->m_btNewLabel->setEnabled(workingData != NULL);
    m_Controls->m_btSaveLabelSet->setEnabled(workingData != NULL);
    m_Controls->m_LabelSetTableWidget->setEnabled(workingData != NULL);
  }
}

void QmitkSegmentationView::ApplyDisplayOptions(mitk::DataNode* node)
{
  if (!node) return;
/*
  bool isBinary(false);
  node->GetPropertyValue("binary", isBinary);

  if (isBinary)
  {
    node->SetProperty( "outline binary", mitk::BoolProperty::New( this->GetPreferences()->GetBool("draw outline", true)) );
    node->SetProperty( "outline width", mitk::FloatProperty::New( 2.0 ) );
    node->SetProperty( "opacity", mitk::FloatProperty::New( this->GetPreferences()->GetBool("draw outline", true) ? 1.0 : 0.3 ) );
    node->SetProperty( "volumerendering", mitk::BoolProperty::New( this->GetPreferences()->GetBool("volume rendering", false) ) );
  }
*/
}

bool QmitkSegmentationView::CheckForSameGeometry(const mitk::DataNode *node1, const mitk::DataNode *node2) const
{
  bool isSameGeometry(true);

  mitk::Image* image1 = dynamic_cast<mitk::Image*>(node1->GetData());
  mitk::Image* image2 = dynamic_cast<mitk::Image*>(node2->GetData());
  if (image1 && image2)
  {
    mitk::Geometry3D* geo1 = image1->GetGeometry();
    mitk::Geometry3D* geo2 = image2->GetGeometry();

    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetOrigin(), geo2->GetOrigin());
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(0), geo2->GetExtent(0));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(1), geo2->GetExtent(1));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(2), geo2->GetExtent(2));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetSpacing(), geo2->GetSpacing());
    isSameGeometry = isSameGeometry && mitk::MatrixEqualElementWise(geo1->GetIndexToWorldTransform()->GetMatrix(), geo2->GetIndexToWorldTransform()->GetMatrix());

    return isSameGeometry;
  }
  else
  {
    return false;
  }
}

void QmitkSegmentationView::UpdateWarningLabel(QString text)
{
  if (text.size() == 0)
    m_Controls->lblSegmentationWarnings->hide();
  else
    m_Controls->lblSegmentationWarnings->show();
  m_Controls->lblSegmentationWarnings->setText(text);
}

void QmitkSegmentationView::CreateQtPartControl(QWidget* parent)
{
  MITK_INFO << "CreateQtPartControl 1";
  // setup the basic GUI of this view
  m_Parent = parent;

  m_Controls = new Ui::QmitkSegmentationControls;
  m_Controls->setupUi(parent);

  m_Controls->patImageSelector->SetDataStorage(this->GetDefaultDataStorage());
  m_Controls->patImageSelector->SetPredicate(m_IsNotLabelSetImagePredicate);

  this->UpdateWarningLabel("Please load an image");

  if( m_Controls->patImageSelector->GetSelectedNode().IsNotNull() )
    this->UpdateWarningLabel("Create a segmentation");

  m_Controls->segImageSelector->SetDataStorage(this->GetDefaultDataStorage());
  m_Controls->segImageSelector->SetPredicate(m_IsOfTypeLabelSetImagePredicate);
  m_Controls->segImageSelector->SetAutoSelectNewItems(true);
  if( m_Controls->segImageSelector->GetSelectedNode().IsNotNull() )
    this->UpdateWarningLabel("");

  mitk::ToolManager* toolManager = m_Controls->m_ManualToolSelectionBox2D->GetToolManager();
  toolManager->SetDataStorage( *(this->GetDefaultDataStorage()) );
  assert ( toolManager );

  m_Controls->m_LabelSetTableWidget->Init();

  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(newLabel()), this, SLOT(OnNewLabel()) );
  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(renameLabel(int)), this, SLOT(OnRenameLabel(int)) );
  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(createSurface(int)), this, SLOT(OnCreateSurface(int)) );
  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(combineAndCreateSurface( const QList<QTableWidgetSelectionRange>& )),
      this, SLOT(OnCombineAndCreateSurface( const QList<QTableWidgetSelectionRange>&)) );

  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(createMask(int)), this, SLOT(OnCreateMask(int)) );
  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(combineAndCreateMask( const QList<QTableWidgetSelectionRange>& )),
      this, SLOT(OnCombineAndCreateMask( const QList<QTableWidgetSelectionRange>&)) );

  m_Controls->m_LabelSearchBox->setAlwaysShowClearIcon(true);
  m_Controls->m_LabelSearchBox->setShowSearchIcon(true);

  QStringList completionList;
  completionList << "";
  m_Completer = new QCompleter(completionList,m_Parent);
  m_Completer->setCaseSensitivity(Qt::CaseInsensitive);
  m_Controls->m_LabelSearchBox->setCompleter(m_Completer);

  connect( m_Controls->m_LabelSearchBox, SIGNAL(returnPressed()), this, SLOT(OnSearchLabel()) );
  connect( m_Controls->m_LabelSetTableWidget, SIGNAL(labelListModified(const QStringList&)), this, SLOT( OnLabelListModified(const QStringList&)) );

  QStringListModel* completeModel = static_cast<QStringListModel*> (m_Completer->model());
  completeModel->setStringList(m_Controls->m_LabelSetTableWidget->GetLabelList());

  //use the same ToolManager instance for our 3D Tools
  m_Controls->m_ManualToolSelectionBox3D->SetToolManager(*toolManager);

  // all part of open source MITK
  m_Controls->m_ManualToolSelectionBox2D->SetGenerateAccelerators(true);
  m_Controls->m_ManualToolSelectionBox2D->SetToolGUIArea( m_Controls->m_ManualToolGUIContainer2D );
  m_Controls->m_ManualToolSelectionBox2D->SetDisplayedToolGroups("Add Subtract Correction Paint Wipe 'Region Growing' Fill Erase 'Live Wire' 'FastMarching2D'");
  m_Controls->m_ManualToolSelectionBox2D->SetLayoutColumns(3);
  m_Controls->m_ManualToolSelectionBox2D->SetEnabledMode( QmitkToolSelectionBox::EnabledWithReferenceAndWorkingDataVisible );
  connect( m_Controls->m_ManualToolSelectionBox2D, SIGNAL(ToolSelected(int)), this, SLOT(OnManualTool2DSelected(int)) );

  //setup 3D Tools
  m_Controls->m_ManualToolSelectionBox3D->SetGenerateAccelerators(true);
  m_Controls->m_ManualToolSelectionBox3D->SetToolGUIArea( m_Controls->m_ManualToolGUIContainer3D );
  //specify tools to be added to 3D Tool area
  m_Controls->m_ManualToolSelectionBox3D->SetDisplayedToolGroups("Threshold 'Two Thresholds' Otsu FastMarching3D RegionGrowing Watershed");
  m_Controls->m_ManualToolSelectionBox3D->SetLayoutColumns(3);
  m_Controls->m_ManualToolSelectionBox3D->SetEnabledMode( QmitkToolSelectionBox::EnabledWithReferenceAndWorkingDataVisible );

  //Hide 3D selection box, show 2D selection box
  m_Controls->m_ManualToolSelectionBox3D->hide();
  m_Controls->m_ManualToolSelectionBox2D->show();

  toolManager->NewNodeObjectsGenerated +=
      mitk::MessageDelegate1<QmitkSegmentationView, mitk::ToolManager::DataVectorType*>( this, &QmitkSegmentationView::NewNodeObjectsGenerated );      // update the list of segmentations

  // create signal/slot connections
  connect( m_Controls->patImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
       this, SLOT( OnPatientComboBoxSelectionChanged( const mitk::DataNode* ) ) );
  connect( m_Controls->segImageSelector, SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
       this, SLOT( OnSegmentationComboBoxSelectionChanged( const mitk::DataNode* ) ) );

  connect( m_Controls->m_btNewLabelSet, SIGNAL(clicked()), this, SLOT(OnNewLabelSet()) );
  connect( m_Controls->m_btSaveLabelSet, SIGNAL(clicked()), this, SLOT(OnSaveLabelSet()) );
  connect( m_Controls->m_btLoadLabelSet, SIGNAL(clicked()), this, SLOT(OnLoadLabelSet()) );
  connect( m_Controls->m_btNewLabel, SIGNAL(clicked()), this, SLOT(OnNewLabel()) );

//  connect( m_Controls->m_pbCreateLabelFromSurface, SIGNAL(clicked()), this, SLOT(CreateLabelFromSurface()) );
//  connect( m_Controls->widgetStack, SIGNAL(currentChanged(int)), this, SLOT(ToolboxStackPageChanged(int)) );


  connect( m_Controls->tabWidgetSegmentationTools, SIGNAL(currentChanged(int)), this, SLOT(OnTabWidgetChanged(int)));

//  connect(m_Controls->MaskSurfaces,  SIGNAL( OnSelectionChanged( const mitk::DataNode* ) ),
//      this, SLOT( OnSurfaceSelectionChanged( ) ) );

  connect(m_Controls->m_SlicesInterpolator, SIGNAL(SignalShowMarkerNodes(bool)), this, SLOT(OnShowMarkerNodes(bool)));
  connect(m_Controls->m_SlicesInterpolator, SIGNAL(Signal3DInterpolationEnabled(bool)), this, SLOT(On3DInterpolationEnabled(bool)));

//  m_Controls->MaskSurfaces->SetDataStorage(this->GetDefaultDataStorage());
//  m_Controls->MaskSurfaces->SetPredicate(mitk::NodePredicateDataType::New("Surface"));

  m_Controls->m_btNewLabelSet->setEnabled(false);
  m_Controls->m_btSaveLabelSet->setEnabled(false);
  m_Controls->m_btLoadLabelSet->setEnabled(false);
  m_Controls->m_btNewLabel->setEnabled(false);

  m_Controls->m_LabelSetTableWidget->setEnabled(false);

  MITK_INFO << "CreateQtPartControl N";
}

void QmitkSegmentationView::OnManualTool2DSelected(int id)
{
    if (id >= 0)
    {
        std::string text = "Active Tool: \"";
        mitk::ToolManager* toolManager = m_Controls->m_ManualToolSelectionBox2D->GetToolManager();
        text += toolManager->GetToolById(id)->GetName();
        text += "\"";
        mitk::StatusBar::GetInstance()->DisplayText(text.c_str());

        mitk::ModuleResource resource = toolManager->GetToolById(id)->GetCursorIconResource();
        this->SetMouseCursor(resource, 0, 0);
    }
    else
    {
        this->ResetMouseCursor();
        mitk::StatusBar::GetInstance()->DisplayText("");
    }
}

void QmitkSegmentationView::ResetMouseCursor()
{
  if ( m_MouseCursorSet )
  {
    mitk::ApplicationCursor::GetInstance()->PopCursor();
    m_MouseCursorSet = false;
  }
}

void QmitkSegmentationView::SetMouseCursor( const mitk::ModuleResource resource, int hotspotX, int hotspotY )
{
  // Remove previously set mouse cursor
  if ( m_MouseCursorSet )
  {
    mitk::ApplicationCursor::GetInstance()->PopCursor();
  }

  mitk::ApplicationCursor::GetInstance()->PushCursor( resource, hotspotX, hotspotY );
  m_MouseCursorSet = true;
}

// ATTENTION some methods for handling the known list of (organ names, colors) are defined in QmitkSegmentationOrganNamesHandling.cpp
