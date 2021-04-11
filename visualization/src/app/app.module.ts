import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { CommonModule } from "@angular/common"
import * as plotlyjs from "plotly.js";
//const PlotlyJS = require('plotly.js');
import { PlotlyModule } from "angular-plotly.js"

PlotlyModule.plotlyjs = plotlyjs
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { HomeComponent } from './home/home.component';
import { GraphHomeComponent } from './graph-home/graph-home.component';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import {MatCardModule} from '@angular/material/card';
import {MatButtonModule} from '@angular/material/button';
import { GraphTfComponent } from './graph-tf/graph-tf.component';
import {MatTreeModule} from '@angular/material/tree';
import { MatIconModule} from '@angular/material/icon';
import { MatCheckboxModule } from '@angular/material/checkbox';
import { FormsModule } from '@angular/forms'
import { FlaskApiService } from './service/flask-api.service';
import { FlexLayoutModule } from "@angular/flex-layout";


@NgModule({
  declarations: [
    AppComponent,
    HomeComponent,
    GraphHomeComponent,
    GraphTfComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    BrowserAnimationsModule,
    MatCardModule,
    MatButtonModule,
    MatTreeModule,
    MatIconModule,
    MatCheckboxModule,
    FormsModule,
    FlexLayoutModule,
    PlotlyModule,
    CommonModule
  ],
  providers: [FlaskApiService],
  bootstrap: [AppComponent]
})
export class AppModule { }
