import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { GraphHomeComponent } from './graph-home/graph-home.component';
import { GraphTfComponent } from './graph-tf/graph-tf.component';
import { HomeComponent } from './home/home.component';

const routes: Routes = [
  {path: "home", component: HomeComponent},
  {path: "graph_home", component: GraphHomeComponent},
  {path: "graph_tf", component: GraphTfComponent},


  {path: "", redirectTo: "home", pathMatch: "full"}
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
