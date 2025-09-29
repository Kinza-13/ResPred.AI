# vf/urls.py
from django.urls import path
from .views import (
    vf_predict_view,
    vf_progress_view,
    vf_progress_status,   # JSON polled by progress.html
    vf_status_api,        # optional alias
    vf_outputs_view,      # results page
    vf_outputs_zip,       # download-all ZIPs
)

app_name = "vf"

urlpatterns = [
    # Upload page (GET shows form, POST starts the run and redirects to progress)
    path("upload/", vf_predict_view, name="vf_upload"),
    path("predict/", vf_predict_view, name="vf_predict"),  # alias if you're already using this

    # Progress screen + polling endpoint
    path("progress/<uuid:run_id>/", vf_progress_view, name="vf_progress"),
    path("progress/<uuid:run_id>/status/", vf_progress_status, name="vf_status"),
    path("progress/<uuid:run_id>/status2/", vf_status_api, name="vf_status2"),  # optional alias

    # Results list + ZIP endpoints
    path("results/", vf_outputs_view, name="vf_outputs"),
    path("results/zip/<str:folder>/", vf_outputs_zip, name="vf_outputs_zip"),
]
