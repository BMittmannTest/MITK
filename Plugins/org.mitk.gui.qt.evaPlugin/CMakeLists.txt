project(org_mitk_gui_qt_evaPlugin)

mitk_create_plugin(
  EXPORT_DIRECTIVE EVAPLUGIN_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt
)
