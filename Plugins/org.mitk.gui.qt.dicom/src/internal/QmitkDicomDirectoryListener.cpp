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
#include "QmitkDicomDirectoryListener.h"

#include <QString>
#include <QFile>
#include <mitkLogMacros.h>
#include <QFileInfoList>

QmitkDicomDirectoryListener::QmitkDicomDirectoryListener()
: m_FileSystemWatcher(new QFileSystemWatcher())
{
    connect(m_FileSystemWatcher,SIGNAL(directoryChanged(const QString&)),this,SLOT(OnDirectoryChanged(const QString&))); 
}

QmitkDicomDirectoryListener::~QmitkDicomDirectoryListener()
{
    delete m_FileSystemWatcher;
}


void QmitkDicomDirectoryListener::OnDirectoryChanged(const QString&)
{   
    SetFilesToImport();
    m_ImportingFiles.append(m_FilesToImport);
    emit SignalAddDicomData(m_FilesToImport);
}

void QmitkDicomDirectoryListener::OnDicomImportFinished(const QStringList& finishedFiles)
{
    RemoveFilesFromDirectoryAndImportingFilesList(finishedFiles);
}

void QmitkDicomDirectoryListener::SetFilesToImport()
{   
    m_FilesToImport.clear();
    QDir listenerDirectory(m_DicomListenerDirectory);
    QFileInfoList entries = listenerDirectory.entryInfoList(QDir::Files);
    if(!entries.isEmpty())
    {
        QFileInfoList::const_iterator file;
        for(file = entries.constBegin(); file != entries.constEnd(); ++file )
        {
            if(!m_ImportingFiles.contains((*file).absoluteFilePath()))
            {
                m_FilesToImport.append((*file).absoluteFilePath());
            }
        }
    } 
}

void QmitkDicomDirectoryListener::RemoveFilesFromDirectoryAndImportingFilesList(const QStringList& files)
{
    QStringListIterator fileToDeleteIterator(files);
    while(fileToDeleteIterator.hasNext())
    {
        QFile file(fileToDeleteIterator.next());
        if(m_ImportingFiles.contains(file.fileName()))
        {
            m_ImportingFiles.removeOne(file.fileName());
            file.remove();
        }
    }
}

void QmitkDicomDirectoryListener::SetDicomListenerDirectory(const QString& directory)
{
    if(isOnlyListenedDirectory(directory))
    {
        QDir listenerDirectory = QDir(directory);
        CreateListenerDirectory(listenerDirectory);
        
        m_DicomListenerDirectory=listenerDirectory.absolutePath();
        m_FileSystemWatcher->addPath(m_DicomListenerDirectory);
        MITK_INFO << m_DicomListenerDirectory.toStdString();
    }
}

const QString& QmitkDicomDirectoryListener::GetDicomListenerDirectory()
{
    return m_DicomListenerDirectory;
}

void QmitkDicomDirectoryListener::CreateListenerDirectory(const QDir& directory)
{   
    if(!directory.exists())
    {
        directory.mkpath(directory.absolutePath());
    }
}

bool QmitkDicomDirectoryListener::isOnlyListenedDirectory(const QString& directory)
{
    bool isOnlyListenedDirectory = false;
    if(m_FileSystemWatcher->directories().count()==0||m_FileSystemWatcher->directories().count()==1)
    {
        if(!m_FileSystemWatcher->directories().contains(directory))
        {
            isOnlyListenedDirectory = true;
        }
    }
    return isOnlyListenedDirectory;
}