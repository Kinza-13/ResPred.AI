# amr/urls.py
from django.urls import path
from django.views.generic import RedirectView
from . import views

app_name = "amr"

urlpatterns = [
    # Canonical upload/start page
    path("upload/", views.amr_predict_view, name="amr_predict"),

    # Back-compat + QoL redirects
    path("", RedirectView.as_view(pattern_name="amr:amr_predict", permanent=False)),
    path("predict/", RedirectView.as_view(pattern_name="amr:amr_predict", permanent=False)),
    path("amr_upload/", RedirectView.as_view(pattern_name="amr:amr_predict", permanent=False), name="amr_upload"),

    # Live progress endpoints
    path("progress/<uuid:run_id>/", views.amr_progress_view, name="amr_progress"),
    path("progress/<uuid:run_id>/status/", views.amr_progress_status, name="amr_progress_status"),

    # Results + ZIP
    path("results/", views.amr_outputs_view, name="amr_outputs"),
    path("results/zip/<str:folder>/", views.amr_outputs_zip, name="amr_outputs_zip"),
]
