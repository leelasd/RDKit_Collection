from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()
import views
urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^Viewer/', include('Viewer.urls',namespace="Viewer")),
    url(r'^OOMMPPAA/', include('OOMMPPAA.urls',namespace="OOMMPPAA")),
    url(r'^admin/', include(admin.site.urls)),
    url( r'^upload/', views.upload, name = 'jfu_upload' ),
    url(r'^run/$', views.run, name='run'),
)+static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
