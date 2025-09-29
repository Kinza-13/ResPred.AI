from django import forms

class UploadFastaForm(forms.Form):
    fasta_file = forms.FileField(label="Upload your FASTA file")
