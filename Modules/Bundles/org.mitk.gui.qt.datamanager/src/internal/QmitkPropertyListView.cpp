#include <cherryISelectionProvider.h>
#include <cherryISelectionService.h>
#include <cherryIWorkbenchWindow.h>
#include <cherryISelectionService.h>
#include "QmitkPropertyListView.h"
#include <QmitkPropertiesTableEditor.h>
#include <mitkDataTreeNodeObject.h>
#include <mbilog.h>
#include <QVBoxLayout>

const std::string QmitkPropertyListView::VIEW_ID = "org.mitk.views.propertylistview";

QmitkPropertyListView::QmitkPropertyListView()
: m_SelectionListener(0)
{
}

QmitkPropertyListView::~QmitkPropertyListView()
{
  cherry::ISelectionService* s = GetSite()->GetWorkbenchWindow()->GetSelectionService();
  if(s)
    s->RemoveSelectionListener(m_SelectionListener);

}

void QmitkPropertyListView::CreateQtPartControl( QWidget* parent )
{
  m_NodePropertiesTableEditor = new QmitkPropertiesTableEditor;

  QVBoxLayout* layout = new QVBoxLayout;
  layout->setContentsMargins(0, 1, 0, 0);
  layout->addWidget(m_NodePropertiesTableEditor);

  parent->setLayout(layout);

  m_SelectionListener = new cherry::SelectionChangedAdapter<QmitkPropertyListView>
    (this, &QmitkPropertyListView::SelectionChanged);
  cherry::ISelectionService* s = GetSite()->GetWorkbenchWindow()->GetSelectionService();
  if(s)
    s->AddSelectionListener(m_SelectionListener);

  cherry::ISelection::ConstPointer selection( this->GetSite()->GetWorkbenchWindow()->GetSelectionService()->GetSelection("org.mitk.views.datamanager"));
  if(selection.IsNotNull())
    this->SelectionChanged(cherry::IWorkbenchPart::Pointer(0), selection);

}

void QmitkPropertyListView::SelectionChanged( cherry::IWorkbenchPart::Pointer, cherry::ISelection::ConstPointer selection )
{
  mitk::DataTreeNodeSelection::ConstPointer _DataTreeNodeSelection 
    = selection.Cast<const mitk::DataTreeNodeSelection>();

  if(_DataTreeNodeSelection.IsNotNull())
  {
    std::vector<mitk::DataTreeNode*> selectedNodes;
    mitk::DataTreeNodeObject* _DataTreeNodeObject = 0;

    for(mitk::DataTreeNodeSelection::iterator it = _DataTreeNodeSelection->Begin();
      it != _DataTreeNodeSelection->End(); ++it)
    {
      _DataTreeNodeObject = dynamic_cast<mitk::DataTreeNodeObject*>((*it).GetPointer());
      if(_DataTreeNodeObject)
        selectedNodes.push_back( _DataTreeNodeObject->GetDataTreeNode() );
    }

    if(selectedNodes.size() > 0)
    {
      m_NodePropertiesTableEditor->SetPropertyList(selectedNodes.back()->GetPropertyList());
    }
  }
}
